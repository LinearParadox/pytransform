import os
import sys
from setuptools import setup

setup(
    name="pytransform",  
    version="0.1",
    description="Additions to the automation of Single-Cell Analysis: A Python package with additional features for automated quality control in single-cell analysis",
    author="Alexei Martinskovsky, Drew Hulsy, Brandon Nguyen",  
    author_email="alexei.martsinkovskiy@gmail.com, drewhulsy@gmail.com, bhn002@.ucsd.edu",  
    url="https://github.com/LinearParadox/pytransform",
    package_dir={"":"pytransform"},
    packages=["pytransform", "pytransform/glm"],
    install_requires=[
        'numpy>=1.20.0',
        'scipy>=1.6.0',
        'scanpy>=1.7.2',
        'anndata>=0.7.6',
        'pandas>=1.2.3',
        'statsmodels>=0.12.0',
        'KDEpy>=1.1.3'
    ]
)

