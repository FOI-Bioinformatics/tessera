import setuptools

setuptools.setup(
    name='recomfi',
    version='1.0',
    author='jaclew',
    description='no_description',
    packages=['recomfi'],
    entry_points= { 'console_scripts': ['recomfi=recomfi:main'] },
    scripts=['recomfi/recomfi.py','recomfi/msa.py','recomfi/recomb.py']
)
