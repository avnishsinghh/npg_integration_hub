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

"""Retrieval of final QC outcomes for loading to the warehouse.

A Python port of the Perl npg_warehouse::loader::illumina::fqc module
(npg_ml_warehouse repository). It retrieves the overall, sequencing,
library and user QC outcomes for a set of compositions and combines them
using the same rules as the Perl code:

  - Only final sequencing and library outcomes are used. A non-final
    outcome is treated as if it did not exist.
  - User QC outcomes have no notion of final, all of them are used.
  - A sequencing outcome for a multi-component entity that has no
    outcome of its own is built from the outcomes for its lanes: any
    fail gives a fail, otherwise any undefined gives undefined,
    otherwise a pass.
  - An entity that is not in the QC database inherits the final library
    outcome of a merged entity it is a component of.
  - The overall outcome is the sequencing outcome, turned into the
    library outcome if sequencing passed and the library outcome is
    defined. A defined user QC outcome replaces the overall outcome.

Outcomes are encoded as 1 (pass), 0 (fail) or None (undefined).

User QC outcomes are stored in the uqc_outcome_ent table and are loaded
to the warehouse qc_user column. They record whether the end user will
utilise a product, regardless of the manual QC outcomes. 'Undecided' maps
to None, so the overall outcome then follows the sequencing and library
outcomes.

The module can also be run as a script, which prints the outcomes for
all compositions of a run as JSON:

  python -m npg_integration_hub.qc.final_qc --db-url mysql+pymysql://... --id-run 45678
"""

import argparse
import json
import logging
import os
import re
import sys
from collections.abc import Iterable, Iterator, Mapping, Sequence
from typing import NamedTuple

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
from sqlalchemy import select
from sqlalchemy.orm import Session, joinedload, selectinload

from npg_integration_hub.db.engine import read_only_engine

logger = logging.getLogger(__name__)

COL_NAME_QC = "qc"
COL_NAME_QC_SEQ = f"{COL_NAME_QC}_seq"
COL_NAME_QC_LIB = f"{COL_NAME_QC}_lib"
COL_NAME_QC_USER = f"{COL_NAME_QC}_user"

SEQ_OUTCOME_ENT = "mqc_outcome_ent"
LIB_OUTCOME_ENT = "mqc_library_outcome_ent"
UQC_OUTCOME_ENT = "uqc_outcome_ent"
COMPONENT_LANES = "component_lanes"
IS_SINGLE_LANE = "is_single_lane"

# Maximum number of values in one SQL IN clause.
MAX_NUMBER_DIGESTS = 900

DB_URL_VAR = "NPGQC_DB_URL"

_FINAL = re.compile(r"final\Z")
_ACCEPTED = re.compile(r"\AAccepted")
_REJECTED = re.compile(r"\ARejected")


class Component(NamedTuple):
    """A component of a composition, see npg_tracking::glossary::composition.

    Two components are equal if all their attributes are equal, which is
    the same as their single-component compositions having the same digest.
    """

    id_run: int
    position: int
    tag_index: int | None = None
    subset: str | None = None

    @property
    def lane_key(self) -> str:
        """The rpt key of the lane this component belongs to, ie '22:1'."""
        return f"{self.id_run}:{self.position}"

    @classmethod
    def from_orm(cls, c: SeqComponent) -> "Component":
        return cls(c.id_run, c.position, c.tag_index, c.subset)


def components_of(composition: SeqComposition) -> list[Component]:
    """Returns the components of a composition from the QC database."""
    return [
        Component.from_orm(scc.seq_component)
        for scc in composition.seq_component_composition
    ]


def is_final(short_desc: str) -> bool:
    """True if an mqc outcome description is final, ie 'Accepted final'."""
    return bool(_FINAL.search(short_desc))


def outcome_value(short_desc: str) -> int | None:
    """Returns 1 for an accepted, 0 for a rejected and None for any other outcome."""
    if _ACCEPTED.search(short_desc):
        return 1
    if _REJECTED.search(short_desc):
        return 0
    return None


def _is_single_lane(components: Sequence[Component]) -> int:
    return 1 if len(components) == 1 and components[0].tag_index is None else 0


def _lanes(components: Iterable[Component]) -> list[str]:
    return [c.lane_key for c in components]


def _chunks(values: Sequence, size: int = MAX_NUMBER_DIGESTS) -> Iterator[Sequence]:
    for i in range(0, len(values), size):
        yield values[i : i + size]


class FinalQCOutcomes:
    """Retrieves and caches final QC outcomes for a set of compositions.

    The equivalent of npg_warehouse::loader::illumina::fqc.

    For maximum efficiency, give the constructor all entities you will
    later need QC outcomes for. Outcomes for all of them are retrieved in
    one go and cached for subsequent requests.

    Example:
      fqc = FinalQCOutcomes(session, {
          "digest1": [Component(22, 1)],
          "digest2": [Component(22, 1, 3), Component(22, 2, 3)],
      })
      fqc.retrieve()
      outcomes = fqc.retrieve_outcomes("digest1")
      outcome = fqc.retrieve_seq_outcome("22:1")

    Giving an empty mapping to the constructor is not an error.

    Requesting outcomes for a digest that was not given to the constructor
    is not an error either. All outcomes are returned as undefined, even
    if final outcomes for this entity exist in the database.
    """

    def __init__(
        self,
        session: Session,
        digests: Mapping[str, Sequence[Component | tuple]],
    ):
        """
        Args:
          session :
            A session for the NPG QC database.
          digests :
            Composition digests mapped to the components of these
            compositions. The components can be Component objects or
            (id_run, position, tag_index, subset) tuples, trailing values
            are optional. This mapping replaces the composition objects
            given to the Perl module. The components are needed for
            entities that are not in the QC database.
        """
        self.session = session
        self.digests: dict[str, list[Component]] = {
            d: [Component(*c) for c in comps] for d, comps in digests.items()
        }
        self._outcomes: dict[str, dict] | None = None
        self._lane_seq_outcomes_cache: dict[str, int | None] | None = None

    def retrieve(self) -> None:
        """Eagerly retrieves and caches the outcomes for all digests.

        Optional; outcomes are otherwise retrieved on first use. Call it
        to control when the database query happens. Repeated calls use
        the cache and do not query the database again.
        """
        self._get_outcomes()

    def retrieve_outcomes(self, digest: str) -> dict[str, int | None]:
        """Returns the warehouse QC values for an entity.

        The keys are qc, qc_seq and qc_lib, and also qc_user if a user
        QC outcome is defined. The values are 1, 0 or None.

        A fail on sequencing QC leads to an overall fail. A pass on
        sequencing QC leads to an overall pass unless the library QC
        outcome is a fail, in which case the overall value is a fail. As
        in the Perl code, an undefined library outcome does not make the
        overall value undefined, since legacy data relies on this.

        A defined user QC outcome replaces the overall outcome.
        """
        if not digest:
            raise ValueError("Composition digest is required")
        if not isinstance(digest, str):
            raise TypeError("Digest should be a string")

        h: dict[str, int | None] = {
            COL_NAME_QC: None,
            COL_NAME_QC_SEQ: None,
            COL_NAME_QC_LIB: None,
        }
        outcome = self._get_outcomes().get(digest)
        if outcome is None:
            return h

        h[COL_NAME_QC_LIB] = outcome.get(LIB_OUTCOME_ENT)
        h[COL_NAME_QC_SEQ] = outcome.get(SEQ_OUTCOME_ENT)
        if h[COL_NAME_QC_SEQ] is None and not outcome[IS_SINGLE_LANE]:
            h[COL_NAME_QC_SEQ] = self._get_lane_seq_outcome(digest)

        if outcome.get(UQC_OUTCOME_ENT) is not None:
            h[COL_NAME_QC_USER] = outcome[UQC_OUTCOME_ENT]
            h[COL_NAME_QC] = h[COL_NAME_QC_USER]
        else:
            h[COL_NAME_QC] = h[COL_NAME_QC_SEQ]
            if h[COL_NAME_QC] and h[COL_NAME_QC_LIB] is not None:
                h[COL_NAME_QC] = h[COL_NAME_QC_LIB]

        return h

    def retrieve_all_outcomes(self) -> dict[str, dict[str, int | None]]:
        """Returns retrieve_outcomes() for every digest given to the constructor."""
        return {digest: self.retrieve_outcomes(digest) for digest in self.digests}

    def retrieve_seq_outcome(self, rpt_key: str) -> dict[str, int | None]:
        """Returns the final sequencing outcome for a lane, ie '22:1'.

        If the outcome for this lane is not cached, an attempt to retrieve
        and cache it is made.
        """
        if not rpt_key:
            raise ValueError("rpt key is required")
        if not isinstance(rpt_key, str):
            raise TypeError("rpt key should be a string")
        lane_outcomes = self._get_lane_seq_outcomes()
        if rpt_key not in lane_outcomes:
            self._cache_lane_seq_outcome(rpt_key)
        return {COL_NAME_QC_SEQ: lane_outcomes.get(rpt_key)}

    def _get_outcomes(self) -> dict[str, dict]:
        if self._outcomes is None:
            self._outcomes = self._build_outcomes()
        return self._outcomes

    def _get_lane_seq_outcomes(self) -> dict[str, int | None]:
        if self._lane_seq_outcomes_cache is None:
            # As in Perl, a single lane without a final outcome is cached as
            # None, so no further attempt to retrieve its outcome is made.
            self._lane_seq_outcomes_cache = {
                data[COMPONENT_LANES][0]: data.get(SEQ_OUTCOME_ENT)
                for data in self._get_outcomes().values()
                if data[IS_SINGLE_LANE]
            }
        return self._lane_seq_outcomes_cache

    def _fetch_compositions(self) -> list[SeqComposition]:
        """Fetches the compositions with their QC outcomes and components."""
        rows: list[SeqComposition] = []
        for chunk in _chunks(list(self.digests)):
            statement = (
                select(SeqComposition)
                .where(SeqComposition.digest.in_(chunk))
                .options(
                    selectinload(SeqComposition.mqc_outcome_ent).joinedload(
                        MqcOutcomeEnt.mqc_outcome_dict
                    ),
                    selectinload(SeqComposition.mqc_library_outcome_ent).joinedload(
                        MqcLibraryOutcomeEnt.mqc_library_outcome_dict
                    ),
                    selectinload(SeqComposition.uqc_outcome_ent).joinedload(
                        UqcOutcomeEnt.uqc_outcome_dict
                    ),
                    selectinload(SeqComposition.seq_component_composition).joinedload(
                        SeqComponentComposition.seq_component
                    ),
                )
            )
            rows.extend(self.session.scalars(statement).all())
        return rows

    def _build_outcomes(self) -> dict[str, dict]:
        outcomes: dict[str, dict] = {}
        if not self.digests:
            return outcomes

        lib_outcomes_decomposed: dict[Component, int] = {}

        for crow in self._fetch_compositions():
            digest = crow.digest
            entry = outcomes.setdefault(digest, {})

            # The outcome tables have a unique index on id_seq_composition,
            # so each list has at most one element.
            seq_desc = (
                crow.mqc_outcome_ent[0].mqc_outcome_dict.short_desc
                if crow.mqc_outcome_ent
                else None
            )
            lib_desc = (
                crow.mqc_library_outcome_ent[0].mqc_library_outcome_dict.short_desc
                if crow.mqc_library_outcome_ent
                else None
            )
            uqc_desc = (
                crow.uqc_outcome_ent[0].uqc_outcome_dict.short_desc
                if crow.uqc_outcome_ent
                else None
            )

            # Sequencing and library outcomes: final only. A final library
            # outcome can be undefined ('Undecided final').
            for name, desc in (
                (SEQ_OUTCOME_ENT, seq_desc),
                (LIB_OUTCOME_ENT, lib_desc),
            ):
                if desc is None or not is_final(desc):
                    continue
                o = outcome_value(desc)
                entry[name] = o
                if (
                    o is not None
                    and name == LIB_OUTCOME_ENT
                    and len(self.digests[digest]) > 1
                ):
                    self._cache_lib_outcome(lib_outcomes_decomposed, digest, o)

            # User QC outcome: no notion of final.
            if uqc_desc is not None:
                entry[UQC_OUTCOME_ENT] = outcome_value(uqc_desc)

            components = components_of(crow)
            entry[IS_SINGLE_LANE] = _is_single_lane(components)
            entry[COMPONENT_LANES] = _lanes(components)

        for digest, components in self.digests.items():
            if digest in outcomes:
                continue
            # Either this entity has not been qc-ed or it will never be.
            # For example, for a Standard workflow of a NovaSeq run, library
            # QC is performed on merged entities, individual plexes are never
            # qc-ed. If this entity is a component of a merged entity with a
            # final library outcome, that outcome is used.
            entry = outcomes[digest] = {
                IS_SINGLE_LANE: _is_single_lane(components),
                COMPONENT_LANES: _lanes(components),
            }
            if len(components) == 1 and components[0] in lib_outcomes_decomposed:
                entry[LIB_OUTCOME_ENT] = lib_outcomes_decomposed[components[0]]

        return outcomes

    def _cache_lib_outcome(
        self, lib_outcomes_decomposed: dict[Component, int], digest: str, o: int
    ) -> None:
        """Caches the library outcome of a multi-component composition
        against each of its components.
        """
        for component in self.digests[digest]:
            if component not in lib_outcomes_decomposed:
                lib_outcomes_decomposed[component] = o
            elif lib_outcomes_decomposed[component] != o:
                # Unlikely, but just in case... Ambiguous results are fails.
                logger.warning("Conflicting inferred outcomes for %s", component)
                lib_outcomes_decomposed[component] = 0

    def _cache_lane_seq_outcome(self, lane_key: str) -> bool:
        """Caches a final sequencing outcome for a lane if one exists.

        As in Perl, the lane is found through the id_run and position
        columns of the outcome table.
        """
        id_run, position = (int(v) for v in lane_key.split(":"))
        statement = (
            select(MqcOutcomeEnt)
            .options(joinedload(MqcOutcomeEnt.mqc_outcome_dict))
            .where(MqcOutcomeEnt.id_run == id_run, MqcOutcomeEnt.position == position)
            .order_by(MqcOutcomeEnt.id_mqc_outcome_ent)
            .limit(1)
        )
        orow = self.session.scalars(statement).first()
        if orow is not None and is_final(orow.mqc_outcome_dict.short_desc):
            self._get_lane_seq_outcomes()[lane_key] = (
                1 if _ACCEPTED.search(orow.mqc_outcome_dict.short_desc) else 0
            )
            return True
        return False

    def _get_lane_seq_outcome(self, digest: str) -> int | None:
        """Combines the sequencing outcomes of the lanes of an entity."""
        lane_outcomes = self._get_lane_seq_outcomes()
        undef_outcome_present = False
        for lane_key in self._get_outcomes()[digest][COMPONENT_LANES]:
            if lane_key not in lane_outcomes:
                self._cache_lane_seq_outcome(lane_key)
            o = lane_outcomes.get(lane_key)
            if o is None:
                undef_outcome_present = True
            elif o == 0:
                return 0
        return None if undef_outcome_present else 1


def compositions_for_run(
    session: Session,
    id_run: int,
    position: int | None = None,
    include_subsets: bool = False,
) -> dict[str, list[Component]]:
    """Builds the digests argument of FinalQCOutcomes for the compositions
    in the QC database that contain at least one component of a run.

    By default, compositions of sequence subsets (phix, human, yhuman,
    xahuman, etc.) are left out. The Perl warehouse loader never passes
    them to the fqc module: npg_warehouse::loader::illumina::autoqc drops
    them or folds them into the entity without a subset.

    Args:
      session :
        A session for the NPG QC database.
      id_run :
        NPG tracking run ID.
      position :
        Optional lane number, restricts the compositions to those with
        a component from this lane.
      include_subsets :
        If True, compositions of sequence subsets are included.
    """
    statement = (
        select(SeqComposition)
        .join(SeqComposition.seq_component_composition)
        .join(SeqComponentComposition.seq_component)
        .where(SeqComponent.id_run == id_run)
        .options(
            selectinload(SeqComposition.seq_component_composition).joinedload(
                SeqComponentComposition.seq_component
            )
        )
        .order_by(SeqComposition.id_seq_composition)
        .distinct()
    )
    if position is not None:
        statement = statement.where(SeqComponent.position == position)
    compositions = session.scalars(statement).all()
    digests = {comp.digest: components_of(comp) for comp in compositions}
    if not include_subsets:
        digests = {
            digest: components
            for digest, components in digests.items()
            if all(c.subset is None for c in components)
        }
    return digests


def final_qc_for_run(
    session: Session,
    id_run: int,
    position: int | None = None,
    include_subsets: bool = False,
    include_indexed_lanes: bool = False,
) -> dict[str, dict]:
    """Returns final QC outcomes for all compositions of a run.

    Each digest is mapped to a dictionary with the components of the
    composition and the outcomes from FinalQCOutcomes.retrieve_outcomes.
    Compositions of sequence subsets are left out unless include_subsets
    is True, see compositions_for_run.

    As in npg_warehouse::loader::illumina::product, entities for whole
    indexed lanes are not end products and are left out unless
    include_indexed_lanes is True. A lane is indexed if either
      - any single plex of this lane is in the QC database, or
      - the lane's tag_metrics or tag_decode_stats autoqc result gives
        a non-zero tags_decode_percent, see tags_decode_percent.
    Indexed lanes are still given to FinalQCOutcomes, as in Perl, so
    their sequencing outcomes are cached.
    """
    digests = compositions_for_run(session, id_run, position, include_subsets)
    fqc = FinalQCOutcomes(session, digests)
    fqc.retrieve()

    products = digests
    if not include_indexed_lanes:
        indexed_lanes = _indexed_lanes(digests.values())
        lane_digests = [d for d, comps in digests.items() if _is_single_lane(comps)]
        decoded_lanes = lanes_with_decoded_tags(session, lane_digests)
        products = {
            digest: components
            for digest, components in digests.items()
            if not (
                _is_single_lane(components)
                and (
                    digest in decoded_lanes
                    or components[0].lane_key in indexed_lanes
                )
            )
        }

    return {
        digest: {
            "components": [c._asdict() for c in components],
            **fqc.retrieve_outcomes(digest),
        }
        for digest, components in products.items()
    }


def _indexed_lanes(compositions: Iterable[Sequence[Component]]) -> set[str]:
    """Returns the lane keys of lanes that have a single-plex composition,
    see _indexed_lanes_hash in npg_warehouse::loader::illumina::product.
    """
    return {
        components[0].lane_key
        for components in compositions
        if len(components) == 1 and components[0].tag_index is not None
    }


def _json_hash(value: str | None) -> dict | None:
    """Deserialises a hash column of an autoqc table, stored as JSON."""
    return None if value is None else json.loads(value)


def _num(value) -> float:
    """Converts a count to a number. Counts are often stored as strings,
    ie "6549030". As in Perl, undefined and non-numeric values are 0.
    """
    if value is None or value == "":
        return 0
    try:
        return float(value)
    except (TypeError, ValueError):
        return 0


def _sum_counts(counts: Mapping | None) -> float:
    """Sums the values, as total_reads_count in
    npg_qc::autoqc::role::tag_metrics.
    """
    return sum(_num(v) for v in (counts or {}).values())


def tag_metrics_decode_percent(row: TagMetrics) -> float | None:
    """Returns tags_decode_percent for a tag_metrics result, as
    _tag_metrics in npg_warehouse::loader::illumina::autoqc.

    It is the percentage of PF reads with a perfect match or one mismatch
    to a tag (perfect_matches_percent + one_mismatch_percent in
    npg_qc::autoqc::role::tag_metrics). None if matches_pf_percent is
    not defined, 0 if there are no PF reads.
    """
    if row.matches_pf_percent is None:
        return None
    total = _sum_counts(_json_hash(row.reads_pf_count))
    if not total:
        return 0
    matched = _sum_counts(_json_hash(row.perfect_matches_pf_count)) + _sum_counts(
        _json_hash(row.one_mismatch_matches_pf_count)
    )
    return matched * 100 / total


def tag_decode_stats_decode_percent(row: TagDecodeStats) -> float | None:
    """Returns tags_decode_percent for a tag_decode_stats result, as
    _tag_decode_stats in npg_warehouse::loader::illumina::autoqc.

    It is 100 minus the percentage of good reads without a tag match
    (decoding_perc_good in npg_qc::autoqc::role::tag_decode_stats). As in
    Perl, the percentage of reads without a match is rounded to two
    decimal places first. None if there are no counts or no count for
    reads without a match.
    """
    errors_good = _json_hash(row.errors_good)
    if errors_good is None:
        return None
    total = _sum_counts(errors_good)
    if total == 0 or "no_match" not in errors_good:
        return None
    no_match = _num(errors_good["no_match"])
    return 100 - float(f"{no_match * 100 / total:.2f}")


def tags_decode_percent(
    tag_metrics: TagMetrics | None, tag_decode_stats: TagDecodeStats | None
) -> float | None:
    """Returns tags_decode_percent for a lane from its autoqc results.

    As in npg_warehouse::loader::illumina::autoqc, a tag_decode_stats
    result is ignored if the lane also has a tag_metrics result.
    """
    if tag_metrics is not None:
        return tag_metrics_decode_percent(tag_metrics)
    if tag_decode_stats is not None:
        return tag_decode_stats_decode_percent(tag_decode_stats)
    return None


def lanes_with_decoded_tags(session: Session, digests: Sequence[str]) -> set[str]:
    """Returns the digests, out of the given lane digests, of lanes with a
    non-zero tags_decode_percent.
    """
    decoded: set[str] = set()
    for chunk in _chunks(list(digests)):
        statement = (
            select(SeqComposition)
            .where(SeqComposition.digest.in_(chunk))
            .options(
                selectinload(SeqComposition.tag_metrics),
                selectinload(SeqComposition.tag_decode_stats),
            )
        )
        for comp in session.scalars(statement).all():
            # Both tables have a unique index on id_seq_composition.
            percent = tags_decode_percent(
                comp.tag_metrics[0] if comp.tag_metrics else None,
                comp.tag_decode_stats[0] if comp.tag_decode_stats else None,
            )
            if percent:
                decoded.add(comp.digest)
    return decoded


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Print final QC outcomes for the compositions of a run as JSON.",
    )
    parser.add_argument(
        "--db-url",
        default=os.environ.get(DB_URL_VAR),
        help="SQLAlchemy URL for the NPG QC database, "
        f"defaults to the value of {DB_URL_VAR}",
    )
    parser.add_argument("--id-run", type=int, required=True, help="NPG run ID")
    parser.add_argument("--position", type=int, help="Optional lane number")
    parser.add_argument(
        "--include-subsets",
        action="store_true",
        help="Include compositions of sequence subsets, ie phix, human",
    )
    parser.add_argument(
        "--include-indexed-lanes",
        action="store_true",
        help="Include whole-lane entities of indexed lanes",
    )
    args = parser.parse_args(argv)

    if not args.db_url:
        parser.error(f"--db-url is required if {DB_URL_VAR} is not set")

    logging.basicConfig(level=logging.WARNING)
    engine = read_only_engine(args.db_url)
    try:
        with Session(engine) as session:
            result = final_qc_for_run(
                session,
                args.id_run,
                args.position,
                include_subsets=args.include_subsets,
                include_indexed_lanes=args.include_indexed_lanes,
            )
    finally:
        engine.dispose()

    json.dump(result, sys.stdout, indent=2)
    sys.stdout.write("\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
