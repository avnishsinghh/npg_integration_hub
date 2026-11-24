# -*- coding: utf-8 -*-
#
# Copyright © 2026 Genome Research Ltd. All rights reserved.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""Tests of final QC outcome retrieval against a MySQL database, see conftest.py.

Some tests change the database step by step and create a new
FinalQCOutcomes object after each change.
"""

import datetime
import hashlib
import json
import logging

import pytest
from npgqc.db.schema import (
    MqcLibraryOutcomeEnt,
    MqcOutcomeEnt,
    SeqComponent,
    SeqComponentComposition,
    SeqComposition,
    TagDecodeStats,
    TagMetrics,
    UqcOutcomeEnt,
)
from pytest import mark as m

from npg_integration_hub.qc.final_qc import (
    Component,
    FinalQCOutcomes,
    compositions_for_run,
    final_qc_for_run,
    main,
    tag_decode_stats_decode_percent,
    tag_metrics_decode_percent,
    tags_decode_percent,
)

ID_RUN = 45678
NEW_ID_RUN = 45680

LANE_DIGEST = "a" * 64  # 45678:1 + 45678:2, Accepted final seq
PLEX_DIGEST = "b" * 64  # 45678:1:3, Accepted final lib, Rejected uqc
PLEX_L1_T5_DIGEST = "c" * 64  # 45678:1:5, Rejected final lib
PLEX_L2_T4_DIGEST = "d" * 64  # 45678:2:4, Accepted preliminary lib, Accepted uqc
OTHER_RUN_DIGEST = "e" * 64  # 45679:1:1, Accepted final lib, Accepted uqc

FIXTURE_DIGESTS = {
    LANE_DIGEST: [Component(ID_RUN, 1), Component(ID_RUN, 2)],
    PLEX_DIGEST: [Component(ID_RUN, 1, 3)],
    PLEX_L1_T5_DIGEST: [Component(ID_RUN, 1, 5)],
    PLEX_L2_T4_DIGEST: [Component(ID_RUN, 2, 4)],
    OTHER_RUN_DIGEST: [Component(45679, 1, 1)],
}

# Outcome dictionary ids, see tests/data/fixtures/*_outcome_dict.yml. The
# ids are the same for sequencing and library outcomes.
ACCEPTED_PRELIMINARY = 1
ACCEPTED_FINAL = 3
REJECTED_FINAL = 4
UNDECIDED_FINAL = 6  # library outcomes only
UQC_ACCEPTED = 1
UQC_REJECTED = 2
UQC_UNDECIDED = 3

CONFLICT_MESSAGE = "Conflicting inferred outcomes for"

# Hash columns of autoqc tables are stored as JSON, often with string counts.
TM_DECODED = {
    "reads_pf_count": json.dumps({"1": "50", "2": "40", "0": "30"}),
    "perfect_matches_pf_count": json.dumps({"1": "48", "2": "42", "0": "0"}),
    "one_mismatch_matches_pf_count": json.dumps({"1": "2", "2": "4", "0": None}),
    "matches_pf_percent": json.dumps({"1": "0.4", "2": "0.35"}),
}
TM_ZERO = {
    **TM_DECODED,
    "perfect_matches_pf_count": json.dumps({"1": "0", "2": "0"}),
    "one_mismatch_matches_pf_count": json.dumps({"1": "0", "2": "0"}),
}
TDS_DECODED = {
    "errors_good": json.dumps({"1": "1505", "0": "275562", "no_match": "2078"}),
}


def _outcomes(qc, qc_seq, qc_lib, **kwargs):
    return {"qc": qc, "qc_seq": qc_seq, "qc_lib": qc_lib, **kwargs}


def _undefined():
    return _outcomes(None, None, None)


########################################################################
# Database helpers
########################################################################
def _freeze(component: Component) -> str:
    """Canonical JSON of a component, undefined attributes left out."""
    return json.dumps(
        {k: v for k, v in component._asdict().items() if v is not None},
        sort_keys=True,
        separators=(",", ":"),
    )


def _digest(*components: Component) -> str:
    """The digest of a composition, as stored in seq_composition.digest.

    It is the SHA256 hex digest of the canonical JSON of the composition,
    with its components sorted by their canonical JSON.
    """
    frozen = sorted(_freeze(c) for c in components)
    data = '{"components":[' + ",".join(frozen) + "]}"
    return hashlib.sha256(data.encode()).hexdigest()


def _add_composition(session, digest, components):
    """Adds a composition with its components to the session."""
    size = len(components)
    compos = SeqComposition(digest=digest, size=size)
    session.add(compos)
    session.flush()
    for c in components:
        comp = (
            session.query(SeqComponent)
            .filter_by(
                id_run=c.id_run,
                position=c.position,
                tag_index=c.tag_index,
                subset=c.subset,
            )
            .one_or_none()
        )
        if comp is None:
            comp = SeqComponent(
                id_run=c.id_run,
                position=c.position,
                tag_index=c.tag_index,
                subset=c.subset,
                digest=hashlib.sha256(_freeze(c).encode()).hexdigest(),
            )
            session.add(comp)
            session.flush()
        session.add(
            SeqComponentComposition(
                id_seq_component=comp.id_seq_component,
                id_seq_composition=compos.id_seq_composition,
                size=size,
            )
        )
    session.flush()
    return compos


def _find_or_save_composition(session, *components: Component) -> SeqComposition:
    """Returns the composition of the components, saving it if needed."""
    digest = _digest(*components)
    compos = session.query(SeqComposition).filter_by(digest=digest).one_or_none()
    return compos or _add_composition(session, digest, components)


def _sync(session):
    """Writes pending changes and expires loaded objects, so that the next
    FinalQCOutcomes object reads the database afresh.
    """
    session.flush()
    session.expire_all()


def _add_lib_outcome(session, compos, id_outcome) -> MqcLibraryOutcomeEnt:
    row = MqcLibraryOutcomeEnt(
        id_seq_composition=compos.id_seq_composition,
        id_mqc_outcome=id_outcome,
        username="cat",
        modified_by="cat",
    )
    session.add(row)
    _sync(session)
    return row


def _add_seq_outcome(session, compos, id_outcome, id_run, position) -> MqcOutcomeEnt:
    row = MqcOutcomeEnt(
        id_seq_composition=compos.id_seq_composition,
        id_mqc_outcome=id_outcome,
        id_run=id_run,
        position=position,
        username="cat",
        modified_by="cat",
    )
    session.add(row)
    _sync(session)
    return row


def _add_uqc_outcome(session, id_seq_composition, id_outcome) -> UqcOutcomeEnt:
    row = UqcOutcomeEnt(
        id_seq_composition=id_seq_composition,
        id_uqc_outcome=id_outcome,
        rationale="RT#456789",
        username="cat",
        modified_by="cat",
    )
    session.add(row)
    _sync(session)
    return row


def _save_lane_outcome(session, id_run, position, id_outcome) -> MqcOutcomeEnt:
    compos = _find_or_save_composition(session, Component(id_run, position))
    return _add_seq_outcome(session, compos, id_outcome, id_run, position)


def _save_lib_outcome(session, id_outcome, *components) -> MqcLibraryOutcomeEnt:
    compos = _find_or_save_composition(session, *components)
    return _add_lib_outcome(session, compos, id_outcome)


def _update(session, row, **values):
    for name, value in values.items():
        setattr(row, name, value)
    _sync(session)


def _delete(session, row):
    session.delete(row)
    _sync(session)


@pytest.fixture
def run3_lanes(session):
    """Lanes 1 to 6 of run 3 with sequencing outcomes: lane 2 is accepted
    final, lane 5 is rejected final, the rest are accepted preliminary.
    Returns the digests of the lanes mapped to their components.
    """
    digests = {}
    rows = {}
    for p in range(1, 7):
        rows[p] = _save_lane_outcome(session, 3, p, ACCEPTED_PRELIMINARY)
        digests[_digest(Component(3, p))] = [Component(3, p)]
    now = datetime.datetime.now()
    _update(session, rows[2], id_mqc_outcome=ACCEPTED_FINAL, reported=now)
    _update(session, rows[5], id_mqc_outcome=REJECTED_FINAL, reported=now)
    return digests


########################################################################
# Tests
########################################################################
@m.describe("Composition digests used in tests")
class TestDigest:
    @m.it("Do not depend on the order of components")
    def test_digest(self):
        c1 = Component(4, 1, 3)
        c2 = Component(4, 3, 3)
        assert _digest(c1, c2) == _digest(c2, c1)
        assert _digest(Component(3, 1)) == hashlib.sha256(
            b'{"components":[{"id_run":3,"position":1}]}'
        ).hexdigest()


@m.describe("Computing tags_decode_percent")
class TestTagsDecodePercent:
    @m.it("Computes the percentage of matched PF reads from tag_metrics")
    def test_tag_metrics(self):
        # 90 perfect + 6 one mismatch out of 120 PF reads
        assert tag_metrics_decode_percent(TagMetrics(**TM_DECODED)) == 80

    @m.it("Gives 0 for tag_metrics without PF reads")
    def test_tag_metrics_no_reads(self):
        tm = TagMetrics(**{**TM_DECODED, "reads_pf_count": json.dumps({"1": "0"})})
        assert tag_metrics_decode_percent(tm) == 0

    @m.it("Gives None for tag_metrics without matches_pf_percent")
    def test_tag_metrics_undefined(self):
        tm = TagMetrics(**{**TM_DECODED, "matches_pf_percent": None})
        assert tag_metrics_decode_percent(tm) is None

    @m.it("Computes 100 minus the rounded no_match percentage from tag_decode_stats")
    def test_tag_decode_stats(self):
        # Real tag_decode_stats counts: no_match is 2078 of 279145,
        # 0.744...% rounded to 0.74.
        tds = TagDecodeStats(**TDS_DECODED)
        assert tag_decode_stats_decode_percent(tds) == pytest.approx(99.26)

    @m.it("Gives 0 for tag_decode_stats where no_match rounds to 100%")
    def test_tag_decode_stats_all_unmatched(self):
        tds = TagDecodeStats(errors_good=json.dumps({"0": "1", "no_match": "99999"}))
        assert tag_decode_stats_decode_percent(tds) == 0

    @m.it("Gives None for tag_decode_stats without a no_match count")
    def test_tag_decode_stats_undefined(self):
        tds = TagDecodeStats(errors_good=json.dumps({"0": "10"}))
        assert tag_decode_stats_decode_percent(tds) is None
        assert tag_decode_stats_decode_percent(TagDecodeStats()) is None

    @m.it("Prefers tag_metrics over tag_decode_stats")
    def test_precedence(self):
        tm = TagMetrics(**TM_ZERO)
        tds = TagDecodeStats(**TDS_DECODED)
        assert tags_decode_percent(tm, tds) == 0
        assert tags_decode_percent(None, tds) == pytest.approx(99.26)
        assert tags_decode_percent(None, None) is None


@m.describe("Validating arguments")
class TestArguments:
    @m.it("Accepts an empty set of digests")
    def test_empty(self, session):
        fqc = FinalQCOutcomes(session, {})
        fqc.retrieve()
        assert fqc.retrieve_all_outcomes() == {}

    @m.it("Is instantiated with digests")
    def test_init(self, session, run3_lanes):
        assert isinstance(FinalQCOutcomes(session, run3_lanes), FinalQCOutcomes)

    @m.it("Requires a digest")
    def test_no_digest(self, session):
        fqc = FinalQCOutcomes(session, {})
        for digest in ("", None):
            with pytest.raises(ValueError, match="Composition digest is required"):
                fqc.retrieve_outcomes(digest)

    @m.it("Requires a string digest")
    def test_digest_type(self, session):
        fqc = FinalQCOutcomes(session, {})
        with pytest.raises(TypeError):
            fqc.retrieve_outcomes(["a" * 64])

    @m.it("Requires an rpt key")
    def test_no_rpt_key(self, session):
        fqc = FinalQCOutcomes(session, {})
        for rpt_key in ("", None):
            with pytest.raises(ValueError, match="rpt key is required"):
                fqc.retrieve_seq_outcome(rpt_key)

    @m.it("Accepts tuples as components")
    def test_tuples(self, session):
        fqc = FinalQCOutcomes(session, {PLEX_DIGEST: [(ID_RUN, 1, 3)]})
        assert fqc.digests[PLEX_DIGEST] == [Component(ID_RUN, 1, 3)]


@m.describe("Retrieving outcomes for entities in the QC database")
class TestFixtureOutcomes:
    @m.it("Uses a final sequencing outcome of a multi-lane entity")
    def test_lane(self, session):
        fqc = FinalQCOutcomes(session, FIXTURE_DIGESTS)
        assert fqc.retrieve_outcomes(LANE_DIGEST) == _outcomes(1, 1, None)

    @m.it("Lets a user QC outcome override a sequencing and library pass")
    def test_uqc_reject(self, session):
        fqc = FinalQCOutcomes(session, FIXTURE_DIGESTS)
        assert fqc.retrieve_outcomes(PLEX_DIGEST) == _outcomes(0, 1, 1, qc_user=0)

    @m.it("Fails an entity with a failed final library outcome")
    def test_lib_reject(self, session):
        fqc = FinalQCOutcomes(session, FIXTURE_DIGESTS)
        assert fqc.retrieve_outcomes(PLEX_L1_T5_DIGEST) == _outcomes(0, 1, 0)

    @m.it("Ignores a preliminary library outcome")
    def test_lib_preliminary(self, session):
        fqc = FinalQCOutcomes(session, FIXTURE_DIGESTS)
        assert fqc.retrieve_outcomes(PLEX_L2_T4_DIGEST) == _outcomes(
            1, None, None, qc_user=1
        )

    @m.it("Leaves the sequencing outcome undefined for a lane without one")
    def test_no_lane_outcome(self, session):
        fqc = FinalQCOutcomes(session, FIXTURE_DIGESTS)
        assert fqc.retrieve_outcomes(OTHER_RUN_DIGEST) == _outcomes(
            1, None, 1, qc_user=1
        )

    @m.it("Returns undefined outcomes for a digest not given to the constructor")
    def test_unknown_digest(self, session):
        fqc = FinalQCOutcomes(session, {LANE_DIGEST: FIXTURE_DIGESTS[LANE_DIGEST]})
        assert fqc.retrieve_outcomes(PLEX_DIGEST) == _undefined()

    @m.it("Retrieves a sequencing outcome for a lane")
    def test_seq_outcome(self, session):
        fqc = FinalQCOutcomes(session, {})
        assert fqc.retrieve_seq_outcome(f"{ID_RUN}:1") == {"qc_seq": 1}
        assert fqc.retrieve_seq_outcome(f"{ID_RUN}:2") == {"qc_seq": None}


@m.describe("Retrieving outcomes for entities not in the QC database")
class TestInferredOutcomes:
    @m.it("Combines lane outcomes: pass if all lanes passed")
    def test_lanes_pass(self, session):
        digest = "1" * 64
        fqc = FinalQCOutcomes(
            session, {digest: [Component(ID_RUN, 1, 3), Component(ID_RUN, 1, 5)]}
        )
        assert fqc.retrieve_outcomes(digest) == _outcomes(1, 1, None)

    @m.it("Combines lane outcomes: undefined if a lane has no outcome")
    def test_lanes_undef(self, session):
        digest = "1" * 64
        fqc = FinalQCOutcomes(
            session, {digest: [Component(ID_RUN, 1, 3), Component(ID_RUN, 2, 3)]}
        )
        assert fqc.retrieve_outcomes(digest) == _undefined()

    @m.it("Combines lane outcomes: fail if a lane failed")
    def test_lanes_fail(self, session):
        lane = _add_composition(session, "2" * 64, [Component(NEW_ID_RUN, 1)])
        _add_seq_outcome(session, lane, REJECTED_FINAL, NEW_ID_RUN, 1)
        digest = "1" * 64
        fqc = FinalQCOutcomes(
            session,
            {digest: [Component(NEW_ID_RUN, 1, 3), Component(ID_RUN, 2, 3)]},
        )
        assert fqc.retrieve_outcomes(digest) == _outcomes(0, 0, None)
        assert fqc.retrieve_seq_outcome(f"{NEW_ID_RUN}:1") == {"qc_seq": 0}

    @m.it("Gives components the library outcome of a merged entity")
    def test_merged_lib_outcome(self, session):
        plex1 = Component(NEW_ID_RUN, 1, 7)
        plex2 = Component(NEW_ID_RUN, 2, 7)
        merged = _add_composition(session, "3" * 64, [plex1, plex2])
        _add_lib_outcome(session, merged, REJECTED_FINAL)

        fqc = FinalQCOutcomes(
            session,
            {"3" * 64: [plex1, plex2], "4" * 64: [plex1], "5" * 64: [plex2]},
        )
        assert fqc.retrieve_outcomes("3" * 64)["qc_lib"] == 0
        assert fqc.retrieve_outcomes("4" * 64)["qc_lib"] == 0
        assert fqc.retrieve_outcomes("5" * 64)["qc_lib"] == 0

    # Both orders of outcomes, so that the outcome of the merged entity
    # processed last does not give the expected fail by chance.
    @pytest.mark.parametrize(
        "outcome1, outcome2",
        [(ACCEPTED_FINAL, REJECTED_FINAL), (REJECTED_FINAL, ACCEPTED_FINAL)],
    )
    @m.it("Fails a component of merged entities with conflicting outcomes")
    def test_conflicting_lib_outcomes(self, session, caplog, outcome1, outcome2):
        plex1 = Component(NEW_ID_RUN, 1, 7)
        plex2 = Component(NEW_ID_RUN, 2, 7)
        plex3 = Component(NEW_ID_RUN, 3, 7)
        merged1 = _add_composition(session, "3" * 64, [plex1, plex2])
        _add_lib_outcome(session, merged1, outcome1)
        merged2 = _add_composition(session, "6" * 64, [plex2, plex3])
        _add_lib_outcome(session, merged2, outcome2)

        fqc = FinalQCOutcomes(
            session,
            {
                "3" * 64: [plex1, plex2],
                "6" * 64: [plex2, plex3],
                "4" * 64: [plex1],
                "5" * 64: [plex2],
                "7" * 64: [plex3],
            },
        )
        value = {ACCEPTED_FINAL: 1, REJECTED_FINAL: 0}
        assert fqc.retrieve_outcomes("4" * 64)["qc_lib"] == value[outcome1]
        assert fqc.retrieve_outcomes("5" * 64)["qc_lib"] == 0  # conflict
        assert fqc.retrieve_outcomes("7" * 64)["qc_lib"] == value[outcome2]
        assert CONFLICT_MESSAGE in caplog.text


@m.describe("Retrieving outcomes for lanes")
class TestLanes:
    @m.it("Retrieves final sequencing outcomes for lanes")
    def test_lane_seq_outcome(self, session, run3_lanes):
        fqc = FinalQCOutcomes(session, run3_lanes)
        assert fqc.retrieve_seq_outcome("3:1") == {"qc_seq": None}, (
            "undefined for a non-final outcome"
        )
        assert fqc.retrieve_seq_outcome("3:2") == {"qc_seq": 1}, (
            "1 for accepted final"
        )
        assert fqc.retrieve_seq_outcome("3:5") == {"qc_seq": 0}, (
            "0 for rejected final"
        )
        assert fqc.retrieve_seq_outcome("3:8") == {"qc_seq": None}, (
            "undef for no record"
        )

    @m.it("Retrieves all outcomes for lanes")
    def test_lane_outcomes(self, session, run3_lanes):
        fqc = FinalQCOutcomes(session, run3_lanes)
        assert fqc.retrieve_outcomes(_digest(Component(3, 1))) == _undefined(), (
            "non-final seq outcome, all undefined"
        )
        assert fqc.retrieve_outcomes(_digest(Component(3, 8))) == _undefined(), (
            "no record, all undefined"
        )
        assert fqc.retrieve_outcomes(_digest(Component(3, 2))) == _outcomes(
            1, 1, None
        ), "accepted final seq outcome result"
        assert fqc.retrieve_outcomes(_digest(Component(3, 5))) == _outcomes(
            0, 0, None
        ), "rejected final seq outcome result"


@m.describe("Retrieving outcomes for a one-component plex")
class TestOneComponentPlex:
    @m.it("Combines sequencing, library and user QC outcomes")
    def test_one_component_plex(self, session, run3_lanes):
        digests = run3_lanes

        def new():
            return FinalQCOutcomes(session, digests)

        def plex(position, id_run=3):
            c = Component(id_run, position, 1)
            return _digest(c), c

        outcomes = _undefined()

        d, c = plex(8)
        digests[d] = [c]
        assert new().retrieve_outcomes(d) == outcomes, "no record - all undefined"

        _save_lib_outcome(session, ACCEPTED_PRELIMINARY, c)
        assert new().retrieve_outcomes(d) == outcomes, (
            "no record for seq, no final record for lib - all undefined"
        )

        d, c = plex(7)
        digests[d] = [c]
        assert new().retrieve_outcomes(d) == outcomes, "no record - all undefined"

        _save_lib_outcome(session, ACCEPTED_FINAL, c)
        outcomes["qc_lib"] = 1
        assert new().retrieve_outcomes(d) == outcomes, (
            "lib final pass, lane - no record - overall undefined"
        )

        outcomes["qc_lib"] = None
        d, c = plex(1)
        digests[d] = [c]
        assert new().retrieve_outcomes(d) == outcomes, (
            "no lib outcome, non-final seq outcome - all undefined"
        )

        d, c = plex(3)
        digests[d] = [c]
        assert new().retrieve_outcomes(d) == outcomes, (
            "no lib outcome, non-final seq outcome - all undefined"
        )

        o = _save_lib_outcome(session, ACCEPTED_PRELIMINARY, c)
        assert new().retrieve_outcomes(d) == outcomes, (
            "both lib and seq outcomes non-final - all undefined"
        )

        _delete(session, o)
        _save_lib_outcome(session, ACCEPTED_FINAL, c)
        outcomes["qc_lib"] = 1
        fqc = new()
        assert fqc.retrieve_outcomes(d) == outcomes, (
            "lib final, seq outcomes non-final - overall undefined"
        )

        outcomes = _undefined()
        d, c = plex(5)
        assert fqc.retrieve_outcomes(d) == outcomes, (
            "all undefined - digest not cached"
        )

        digests[d] = [c]
        outcomes.update(qc_seq=0, qc=0)
        assert new().retrieve_outcomes(d) == outcomes, (
            "final fail seq outcome, no lib outcome - overall fail"
        )

        o = _save_lib_outcome(session, ACCEPTED_PRELIMINARY, c)
        assert new().retrieve_outcomes(d) == outcomes, (
            "seq outcome final, lib not final - overall fail"
        )

        _delete(session, o)
        o = _save_lib_outcome(session, ACCEPTED_FINAL, c)
        outcomes.update(qc_lib=1, qc=0)
        assert new().retrieve_outcomes(d) == outcomes, (
            "seq and lib outcomes final, but opposite - overall fail"
        )

        _delete(session, o)
        o = _save_lib_outcome(session, REJECTED_FINAL, c)
        outcomes.update(qc_lib=0, qc=0)
        assert new().retrieve_outcomes(d) == outcomes, (
            "seq and lib outcomes final all fail - all fail"
        )

        _delete(session, o)
        o = _save_lib_outcome(session, UNDECIDED_FINAL, c)
        outcomes.update(qc_lib=None, qc=0)
        assert new().retrieve_outcomes(d) == outcomes, (
            "seq and lib outcomes final, lib undef - overall fail"
        )

        d, c = plex(2)
        digests[d] = [c]
        outcomes.update(qc_seq=1, qc=1)
        assert new().retrieve_outcomes(d) == outcomes, (
            "final seq outcome, no lib outcome - overall pass"
        )

        o = _save_lib_outcome(session, ACCEPTED_PRELIMINARY, c)
        assert new().retrieve_outcomes(d) == outcomes, (
            "seq outcome pass final, lib not final - overall pass"
        )

        _delete(session, o)
        o = _save_lib_outcome(session, ACCEPTED_FINAL, c)
        outcomes.update(qc_lib=1, qc=1)
        assert new().retrieve_outcomes(d) == outcomes, (
            "seq and lib outcomes final - all pass"
        )

        _delete(session, o)
        o = _save_lib_outcome(session, REJECTED_FINAL, c)
        outcomes.update(qc_lib=0, qc=0)
        assert new().retrieve_outcomes(d) == outcomes, (
            "seq and lib outcomes final, lib fail - lib and qc fail"
        )

        _delete(session, o)
        o = _save_lib_outcome(session, UNDECIDED_FINAL, c)
        outcomes.update(qc_lib=None, qc=1)
        assert new().retrieve_outcomes(d) == outcomes, (
            "seq and lib outcomes final, lib undef, qc pass"
        )

        u = _add_uqc_outcome(session, o.id_seq_composition, UQC_ACCEPTED)
        outcomes["qc_user"] = 1
        assert new().retrieve_outcomes(d) == outcomes, "as above, uqc pass"

        _update(session, u, id_uqc_outcome=UQC_REJECTED)
        outcomes.update(qc_user=0, qc=0)
        assert new().retrieve_outcomes(d) == outcomes, "as above, uqc fail"

        _update(session, u, id_uqc_outcome=UQC_UNDECIDED)
        del outcomes["qc_user"]
        outcomes["qc"] = 1
        assert new().retrieve_outcomes(d) == outcomes, (
            "as above, uqc outcome is not set"
        )

    @m.it("Mirrors the user QC outcome in absence of other outcomes")
    def test_user_qc_only(self, session):
        c = Component(30333, 5, 1)
        d = _digest(c)
        compos = _find_or_save_composition(session, c)
        u = _add_uqc_outcome(session, compos.id_seq_composition, UQC_ACCEPTED)

        def retrieve():
            return FinalQCOutcomes(session, {d: [c]}).retrieve_outcomes(d)

        assert retrieve() == _outcomes(1, None, None, qc_user=1), (
            "overall qc is pass to mirror user qc in absence of other qc outcomes"
        )

        _update(session, u, id_uqc_outcome=UQC_REJECTED)
        assert retrieve() == _outcomes(0, None, None, qc_user=0), (
            "overall qc is fail to mirror user qc in absence of other qc outcomes"
        )

        _update(session, u, id_uqc_outcome=UQC_UNDECIDED)
        assert retrieve() == _undefined(), "none of the qc outcomes is defined"


@m.describe("Retrieving outcomes for a multi-component plex")
class TestMultiComponentPlex:
    @m.it("Retrieves outcomes for merged plexes and their components")
    def test_multi_component_plex(self, session, caplog):
        id_run = 4
        digests = {}

        def new():
            return FinalQCOutcomes(session, digests)

        lane_rows = []
        for p in range(1, 5):
            lane_rows.append(_save_lane_outcome(session, id_run, p, ACCEPTED_FINAL))
            digests[_digest(Component(id_run, p))] = [Component(id_run, p)]

        # Plexes 1 to 6 merged across lanes 1 to 4. The library outcome of
        # merged plex i has dictionary id i. The components are added as
        # single-component entities that are not in the QC database.
        compositions = []
        compositions_unmerged = []
        for i in range(1, 7):
            components = [Component(id_run, p, i) for p in range(1, 5)]
            _save_lib_outcome(session, i, *components)
            compositions.append(components)
            digests[_digest(*components)] = components
            for c in components:
                digests[_digest(c)] = [c]
                compositions_unmerged.append([c])

        def tag_index(components):
            return components[0].tag_index

        def lib_outcome(i):
            return 1 if i == 3 else (0 if i == 4 else None)

        fqc = new()
        for c in compositions + compositions_unmerged:
            i = tag_index(c)
            expected = _outcomes(0 if i == 4 else 1, 1, lib_outcome(i))
            assert fqc.retrieve_outcomes(_digest(*c)) == expected, (
                f"correct outcome for plex {i} when all lanes pass"
            )

        for lane in lane_rows:
            _update(session, lane, id_mqc_outcome=REJECTED_FINAL)
        fqc = new()
        for c in compositions + compositions_unmerged:
            i = tag_index(c)
            assert fqc.retrieve_outcomes(_digest(*c)) == _outcomes(
                0, 0, lib_outcome(i)
            ), f"correct outcome for plex {i} when all lanes fail"

        _update(session, lane_rows[0], id_mqc_outcome=ACCEPTED_FINAL)
        fqc = new()
        for c in compositions:
            i = tag_index(c)
            assert fqc.retrieve_outcomes(_digest(*c)) == _outcomes(
                0, 0, lib_outcome(i)
            ), f"correct outcome for {i} when some lanes fail, some lanes pass"

        for lane in lane_rows:
            _update(session, lane, id_mqc_outcome=ACCEPTED_FINAL)
        _update(session, lane_rows[1], id_mqc_outcome=ACCEPTED_PRELIMINARY)
        fqc = new()
        for c in compositions + compositions_unmerged:
            i = tag_index(c)
            expected = _outcomes(None, None, lib_outcome(i))
            if len(c) == 1 and c[0].position != 2:
                expected.update(qc_seq=1, qc=0 if i == 4 else 1)
            assert fqc.retrieve_outcomes(_digest(*c)) == expected, (
                f"correct outcome for {i} when one lane undef, others pass"
            )

        c = Component(id_run, 1, 4)
        row = _save_lib_outcome(session, ACCEPTED_FINAL, c)
        d = _digest(c)
        outcomes = _outcomes(1, 1, 1)
        assert new().retrieve_outcomes(d) == outcomes, (
            "existing lib value is not overwritten"
        )

        u = _add_uqc_outcome(session, row.id_seq_composition, UQC_ACCEPTED)
        outcomes["qc_user"] = 1
        assert new().retrieve_outcomes(d) == outcomes, "as above, uqc pass, qc pass"

        _update(session, u, id_uqc_outcome=UQC_REJECTED)
        outcomes.update(qc_user=0, qc=0)
        assert new().retrieve_outcomes(d) == outcomes, (
            "as above, uqc fail, leading to qc fail"
        )

        _update(session, u, id_uqc_outcome=UQC_UNDECIDED)
        del outcomes["qc_user"]
        outcomes["qc"] = 1
        assert new().retrieve_outcomes(d) == outcomes, (
            "as above, uqc outcome is not set"
        )

        # Plex 3 from lanes 1 and 3 merged, a library fail. Plex 3 merged
        # across all lanes is a library pass, so the outcomes inferred for
        # 4:1:3 and 4:3:3 are in conflict.
        components = [Component(id_run, 1, 3), Component(id_run, 3, 3)]
        row = _save_lib_outcome(session, REJECTED_FINAL, *components)
        d = _digest(*components)
        digests[d] = components

        def retrieve_with_warnings():
            caplog.clear()
            with caplog.at_level(logging.WARNING, logger="npg_integration_hub.qc.final_qc"):
                retrieved = new().retrieve_outcomes(d)
            warnings = [r for r in caplog.records if CONFLICT_MESSAGE in r.getMessage()]
            return retrieved, warnings

        retrieved, warnings = retrieve_with_warnings()
        assert len(warnings) == 2, "warning about conflicting outcome values"
        assert retrieved == _outcomes(0, 1, 0), (
            "lib and overall outcomes are a fail when inferred lib outcomes "
            "are in conflict"
        )

        _add_uqc_outcome(session, row.id_seq_composition, UQC_ACCEPTED)
        retrieved, warnings = retrieve_with_warnings()
        assert len(warnings) == 2, "warning about conflicting outcome values"
        assert retrieved == _outcomes(1, 1, 0, qc_user=1), (
            "overall fail overwritten by user qc pass"
        )


@m.describe("Retrieving outcomes for a run")
class TestRun:
    @m.it("Finds the components of every composition of a run")
    def test_compositions_for_run(self, session):
        digests = compositions_for_run(session, ID_RUN)
        assert digests == {
            d: c for d, c in FIXTURE_DIGESTS.items() if d != OTHER_RUN_DIGEST
        }

    @m.it("Leaves out compositions of sequence subsets by default")
    def test_subsets_excluded(self, session):
        _add_composition(session, "8" * 64, [Component(ID_RUN, 1, 3, "phix")])
        _add_composition(session, "9" * 64, [Component(ID_RUN, 1, 3, "human")])
        digests = compositions_for_run(session, ID_RUN)
        assert "8" * 64 not in digests
        assert "9" * 64 not in digests
        assert PLEX_DIGEST in digests

    @m.it("Includes compositions of sequence subsets if asked")
    def test_subsets_included(self, session):
        _add_composition(session, "8" * 64, [Component(ID_RUN, 1, 3, "phix")])
        digests = compositions_for_run(session, ID_RUN, include_subsets=True)
        assert digests["8" * 64] == [Component(ID_RUN, 1, 3, "phix")]

    @m.it("Leaves out whole-lane entities of indexed lanes by default")
    def test_indexed_lanes_excluded(self, session):
        lane1 = "8" * 64  # lane 1 has plexes in the fixtures
        lane3 = "9" * 64  # lane 3 has no plexes
        _add_composition(session, lane1, [Component(ID_RUN, 1)])
        _add_composition(session, lane3, [Component(ID_RUN, 3)])
        result = final_qc_for_run(session, ID_RUN)
        assert lane1 not in result
        assert lane3 in result
        assert PLEX_DIGEST in result
        assert LANE_DIGEST in result  # merged lanes are kept

    @m.it("Leaves out lanes with tags decoded according to tag_metrics")
    def test_tag_metrics_lane_excluded(self, session):
        lane3 = _add_composition(session, "9" * 64, [Component(ID_RUN, 3)])
        session.add(TagMetrics(id_seq_composition=lane3.id_seq_composition, **TM_DECODED))
        lane4 = _add_composition(session, "8" * 64, [Component(ID_RUN, 4)])
        session.add(TagMetrics(id_seq_composition=lane4.id_seq_composition, **TM_ZERO))
        session.flush()
        result = final_qc_for_run(session, ID_RUN)
        assert "9" * 64 not in result
        assert "8" * 64 in result  # zero tags decoded, not an indexed lane

    @m.it("Leaves out lanes with tags decoded according to tag_decode_stats")
    def test_tag_decode_stats_lane_excluded(self, session):
        lane3 = _add_composition(session, "9" * 64, [Component(ID_RUN, 3)])
        session.add(
            TagDecodeStats(id_seq_composition=lane3.id_seq_composition, **TDS_DECODED)
        )
        session.flush()
        assert "9" * 64 not in final_qc_for_run(session, ID_RUN)

    @m.it("Ignores tag_decode_stats if a lane has a tag_metrics result")
    def test_tag_metrics_precedence(self, session):
        lane3 = _add_composition(session, "9" * 64, [Component(ID_RUN, 3)])
        id_compos = lane3.id_seq_composition
        session.add(TagMetrics(id_seq_composition=id_compos, **TM_ZERO))
        session.add(TagDecodeStats(id_seq_composition=id_compos, **TDS_DECODED))
        session.flush()
        assert "9" * 64 in final_qc_for_run(session, ID_RUN)

    @m.it("Includes whole-lane entities of indexed lanes if asked")
    def test_indexed_lanes_included(self, session):
        lane1 = "8" * 64
        _add_composition(session, lane1, [Component(ID_RUN, 1)])
        result = final_qc_for_run(session, ID_RUN, include_indexed_lanes=True)
        assert result[lane1]["qc_seq"] is None

    @m.it("Returns outcomes and components for every composition of a run")
    def test_final_qc_for_run(self, session):
        result = final_qc_for_run(session, ID_RUN, 2)
        assert list(result) == [LANE_DIGEST, PLEX_L2_T4_DIGEST]
        assert result[PLEX_L2_T4_DIGEST] == {
            "components": [
                {"id_run": ID_RUN, "position": 2, "tag_index": 4, "subset": None}
            ],
            **_outcomes(1, None, None, qc_user=1),
        }

    @m.it("Prints outcomes as JSON from the command line")
    def test_main(self, loaded_engine, capsys):
        url = loaded_engine.url.render_as_string(hide_password=False)
        assert main(["--db-url", url, "--id-run", str(ID_RUN)]) == 0
        result = json.loads(capsys.readouterr().out)
        assert result[PLEX_DIGEST]["qc"] == 0
        assert result[PLEX_DIGEST]["qc_user"] == 0
