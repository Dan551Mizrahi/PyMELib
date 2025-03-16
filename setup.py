from setuptools import find_packages, setup

setup(
    name='PyMELib',
    packages=find_packages(include=['PyMELib']),
    version='0.1.0',
    description='First version of the PyMELib (Python Minimal Enumeration Library) library',
    author='Dan S. Mizrahi and Batya Kenig',
    install_requires=['networkx', 'typing', 'frozendict'],
    setup_requires=['pytest-runner'],
    tests_require=['pytest'],
    test_suite='tests',
)