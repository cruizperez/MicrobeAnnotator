SHELL=/bin/bash
PKG_NAME=microbeannotator


sortimport:
	isort src/

test/unit:
	pytest -vv -r s src/tests/unit/*.py src/tests/unit/*/*.py --cov=src/${PKG_NAME} --cov-report term-missing

# test/integration:
# 	pytest -vv -r s src/tests/integration/*.py

test/lint:
	flake8 --config=pyproject.toml src/

test/isort:
	isort src/ --check --diff

test/type:
	mypy src

test/style:
	black --line-length 120 --check src

ci: \
	test/lint \
	test/isort \
	test/type \
	test/style \
	test/unit