SHELL=/bin/bash
ENV_NAME=microbeannotator
PKG_NAME=microbeannotator
SNAKEMAKE_FOLDER=src/${PKG_NAME}/.......

CONDA_ACTIVATE=source $$(conda info --base)/etc/profile.d/conda.sh ; conda activate ; conda activate ${ENV_NAME}

install:
    conda-lock install -n ${ENV_NAME} \
    && ${CONDA_ACTIVATE} \
    && pip install -r requirements.txt \
    && python setup.py install

dev_env:
    conda-lock install -n ${ENV_NAME} \
    && ${CONDA_ACTIVATE} \
    && pip install -r requirements.txt \
    && pip install -r requirements.dev.txt \
    && python setup.py develop

lock:
    rm -f conda-lock.yml && conda-lock -f environment.yml -p linux-64

test-unit:
    pytest -vv -x src/tests/unit/*.py --cov=src/${PKG_NAME} --cov-report term-missing

test-integration:
    pytest -vv -x src/tests/integration/*.py

lint:
    flake8 --config=setup.cfg src

typecheck:
    mypy --config=setup.cfg src

stylecheck:
    black --line-length 99 --check src
    snakefmt --line-length 99 --check ${SNAKEMAKE_FOLDER}

ci:
    lint typecheck stylecheck test-unit test-integration
