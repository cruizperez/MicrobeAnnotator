from setuptools import setup, find_packages

setup(
    name='microbeannotator',
    version='1.0.1',
    description='A user friendly microbe genome annotation tool',
    url='https://github.com/cruizperez/MicrobeAnnotator',
    author='Carlos A. Ruiz Perez',
    author_email='cruizperez3@gatech.edu',
    packages=find_packages(),
    include_package_data=True,
    use_scm_version=True,
    setup_requires=[
        'setuptools_scm'],
    python_requires='>=3.5, <3.9',
    scripts=['microbeannotator.py',
             'microbeannotator_db_builder.py', 
             'independent_scripts/identifier_conversion.py',
             'independent_scripts/ko_mapper.py'],
    entry_points={
        'console_scripts': [
            'microbeannotator=microbeannotator:main',
            'microbeannotator_db_builder=microbeannotator_db_builder:main',
        ],
    },
)
