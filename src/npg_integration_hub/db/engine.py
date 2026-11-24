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

"""Database engines, for use with any of the databases this package reads."""

from typing import Any

from sqlalchemy import create_engine, event
from sqlalchemy.engine import Engine


def read_only_engine(url: str, **kwargs: Any) -> Engine:
    """Returns a MySQL engine whose connections cannot write to the database.

    Every new connection is set to SET SESSION TRANSACTION READ ONLY, so
    MySQL rejects writes to non-temporary tables with error 1792. This is
    a client-side safeguard against accidental writes: code on the same
    connection can switch it off again. Use a database account with only
    SELECT privileges where writes must be impossible.

    Args:
      url :
        SQLAlchemy URL for a MySQL database.
      kwargs :
        Further arguments for sqlalchemy.create_engine. pool_pre_ping is
        True unless given.
    """
    kwargs.setdefault("pool_pre_ping", True)
    engine = create_engine(url, **kwargs)

    @event.listens_for(engine, "connect")
    def _set_read_only(dbapi_connection, connection_record):
        with dbapi_connection.cursor() as cursor:
            cursor.execute("SET SESSION TRANSACTION READ ONLY")

    return engine
