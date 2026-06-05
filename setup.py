import sys
import platform
import setuptools
from setuptools import setup

# --- Hard block non-Linux installs ---
if platform.system() != "Linux":
    sys.exit("ERROR: SomatiX is supported only on Linux systems.")

setup(
    name="somatix",
    version="0.0.1",
    author="Zelin Liu",
    author_email="zlliu95@outlook.com",
    description="A long-read DNA variant caller for somatic SNV using tumor-normal paired samples",
    long_description=open("README.md", encoding="utf-8").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/colinliuzelin/SomatiX",
    license="SomatiX Research Use License",
    packages=setuptools.find_packages(where="source"),
    package_dir={"": "source"},
    package_data={
        "somatix": [
            "bin/allele_counter",
            "bin/allele_counter.cpp",
            "bin/Makefile",
        ],
    },
    data_files=[
        ("bin", ["source/somatix/bin/allele_counter"]),
    ],
    include_package_data=False,
    entry_points={
        "console_scripts": [
            "somatix=somatix.somatix:main",
        ]
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.12",
        "Operating System :: POSIX :: Linux",
        "Environment :: Console",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    # Strictly Python 3.12.x
    python_requires=">=3.12,<3.13",
)
