# npg_integration_hub

This package provides the data transformations for the integration hub. Each
transformation retrieves data from a source system and converts it into the
format expected by downstream consumers of the hub.

## Data QC

It retrieves final QC outcomes for sequencing data from the npg_qc database, which
stores both autoqc results and manual QC outcomes. Database access is provided by
the ORM in [npgqc-python](https://github.com/avnishsinghh/npgqc-python).

To output the final QC outcomes for a run in JSON format:

```
python -m npg_integration_hub.qc.final_qc --db-url mysql+pymysql://... --id-run 45678
```

The database URL may alternatively be supplied through the `NPGQC_DB_URL`
environment variable, in which case the `--db-url` option can be omitted. A
value given with `--db-url` takes precedence over the environment variable.
One of the two must be provided.

```
export NPGQC_DB_URL=mysql+pymysql://...
python -m npg_integration_hub.qc.final_qc --id-run 45678
```

## Development

The project follows the Google Python
[style guide](https://github.com/google/styleguide/blob/gh-pages/pyguide.md)
for code and documentation. Linting is performed with `ruff`.

To install the package with its test dependencies and run the test suite:

```
pip install .[test]
pytest
```

Tests that require a MySQL server are skipped unless the `NPGQC_TEST_DB_URL`
environment variable is set. The database named in the URL is created at the
start of the test run and dropped on completion; it must not refer to a
database containing data that needs to be retained.

```
docker run -d --rm --name npgqc-test-mysql -e MYSQL_ROOT_PASSWORD=test -p 13306:3306 mysql:8.0
NPGQC_TEST_DB_URL=mysql+pymysql://root:test@127.0.0.1:13306/npgqc_test pytest
```
