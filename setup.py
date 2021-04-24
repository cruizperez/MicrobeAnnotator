from setuptools import setup, find_packages
from microbeannotator import version
import os

requirements = [req.strip() for req in open('requirements.txt', 'r').readlines()]

with open(os.path.join(os.path.dirname(__file__), 'README.md')) as readme:
    README = readme.read()

setup(
    name='microbeannotator',
    version=version,

    scripts=['bin/microbeannotator', 'bin/microbeannotator_db_builder'],
    include_package_data=True,

    packages=find_packages(),

    install_requires = requirements,
    python_requires='>=3.6, <=3.7.5',
    author='Carlos A. Ruiz Perez',
    author_email='cruizperez3@gatech.edu',
    description='A user friendly microbe genome annotation tool',
    licence='Artistic License 2.0',
    keywords = ['genome annotation', 'protein', 'comparative genomics'],
    url='https://github.com/cruizperez/MicrobeAnnotator'
)
