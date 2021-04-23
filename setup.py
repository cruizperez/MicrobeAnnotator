from setuptools import setup, find_packages
from microbeannotator import version
setup(
    name='microbeannotator',
    version=version,
    description='A user friendly microbe genome annotation tool',
    url='https://github.com/cruizperez/MicrobeAnnotator',
    author='Carlos A. Ruiz Perez',
    author_email='cruizperez3@gatech.edu',
    keywords = ['genome annotation', 'protein', 'comparative genomics'],
    packages=find_packages(),
    include_package_data=True,
    use_scm_version=False,
    setup_requires=[
        'setuptools_scm'],
    python_requires='>=3.5, <3.9',
    scripts=['bin/microbeannotator',
             'bin/microbeannotator_db_builder'] 
)
