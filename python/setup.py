from setuptools import setup
from distutils.core import setup

setup(
    name='fcst',
    version='0.1.0',
    author='FCST development group',
    author_email='developper@fcst.org',
    packages=['fcst', 'fcst.test'],
    scripts=['bin/runALOP.py','bin/runTransient.py'],
    url='http://pypi.python.org/pypi/pythonFCST/',
    license='LICENSE.txt',
    description='Utilities for fcst.',
    long_description=open('README.txt').read(),
    install_requires=[
        "scipy >= 0.9.0",
        "numpy >= 1.6.1",
    ],
)