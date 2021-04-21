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
    packages=['microbeannotator'],#find_packages(),
    include_package_data=True,
    use_scm_version=False,
    setup_requires=[
        'setuptools_scm'],
    python_requires='>=3.5, <3.9',
    scripts=['microbeannotator/microbeannotator.py',
             'microbeannotator/microbeannotator_db_builder.py', 
             'microbeannotator/pipeline/identifier_conversion.py',
             'microbeannotator/pipeline/ko_mapper.py'],
    entry_points={
        'console_scripts': [
            'microbeannotator = microbeannotator:main',
            'microbeannotator_db_builder = microbeannotator_db_builder:main',
        ],
    },
)
