SHELL=/bin/bash
PKG_NAME=microbeannotator


sortimport:
	isort src/

test/unit:
	pytest -vv -x src/tests/unit/*/*.py --cov=src/${PKG_NAME} --cov-report term-missing

test/integration:
	pytest -vv -x src/tests/integration/*.py

test/lint:
	flake8 --config=setup.cfg src

test/isort:
	isort src/ --check --diff

test/type:
	mypy --config=setup.cfg src

test/style:
	black --line-length 99 --check src

ci: \
	test/unit \
	test/lint \
	test/isort \
	test/type \
	test/style \