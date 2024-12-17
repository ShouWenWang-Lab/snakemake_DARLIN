import os
import sys
from pathlib import Path

from setuptools import find_packages, setup

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__))))

setup(
    name="darlin",
    version="0.0.1",
    python_requires=">=3.6",
    packages=find_packages(),  # this is better than packages=["cospar"], which only include the top level files
    author="Shou-Wen Wang",
    author_email="wangshouwen@westlake.edu.cn",
    description="DARLIN snakemake pipeline",
    long_description=Path("README.md").read_text("utf-8"),
    license="BSD",
)
