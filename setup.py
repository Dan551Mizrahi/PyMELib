from setuptools import find_packages, setup

setup(
    name='PyMELib',
    packages=find_packages(include=['PyMELib', 'PyMELib.utils']),
    version='0.6.0',
    description='First version of the PyMELib (Python Minimal Enumeration Library) library',
    author='Dan S. Mizrahi and Batya Kenig',
    author_email="danmizrahithemiz@gmail.com",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    install_requires=['networkx', 'typing', 'frozendict', 'matplotlib', 'EoN', 'plotly'],
    extras_require={
        'test': ['pytest',],
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.8",
)