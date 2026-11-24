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

"""Tests of database engines against a MySQL database, see conftest.py."""

import pytest
from npgqc.db.schema import SeqComposition
from pytest import mark as m
from sqlalchemy import select, update
from sqlalchemy.exc import OperationalError
from sqlalchemy.orm import Session

from npg_integration_hub.db.engine import read_only_engine


@pytest.fixture
def ro_engine(loaded_engine):
    engine = read_only_engine(loaded_engine.url.render_as_string(hide_password=False))
    yield engine
    engine.dispose()


@m.describe("Read-only engine")
class TestReadOnlyEngine:
    @m.it("Reads from the database")
    def test_read(self, ro_engine):
        with Session(ro_engine) as session:
            digests = session.scalars(select(SeqComposition.digest)).all()
        assert "a" * 64 in digests

    @m.it("Rejects writes to the database")
    def test_write(self, ro_engine):
        with Session(ro_engine) as session:
            with pytest.raises(OperationalError, match="1792"):
                session.execute(
                    update(SeqComposition)
                    .where(SeqComposition.id_seq_composition == 1)
                    .values(size=2)
                )

    @m.it("Keeps connections read-only after the pool reconnects")
    def test_reconnect(self, ro_engine):
        with Session(ro_engine) as session:
            session.scalars(select(SeqComposition.digest)).all()
        ro_engine.dispose()
        with ro_engine.connect() as conn:
            assert conn.exec_driver_sql("SELECT @@transaction_read_only").scalar() == 1
