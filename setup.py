import os
import sys
import platform
import setuptools
from setuptools import setup

# --- Hard block non-Linux installs ---
if platform.system() != "Linux":
    sys.exit("ERROR: SomatiX is supported only on Linux systems.")

def read_requirements(filename):
    if os.path.isfile(filename):
        with open(filename) as f:
            return [line.strip() for line in f if line.strip() and not line.startswith("#")]
    return []


# Detect local allele_counter binary and install it as a script if present
def find_allele_counter():
    candidates = [
        os.path.join("bin", "allele_counter"),
        os.path.join("source", "somatix", "bin", "allele_counter"),
    ]
    found = []
    for p in candidates:
        if os.path.isfile(p) and os.access(p, os.X_OK):
            found.append(p)
    return found

_scripts = find_allele_counter()


setup(
    name="somatix",
    version="0.0.1",
    author="Zelin Liu",
    author_email="liuz6@chop.edu",
    description="A long-read DNA variant caller for somatic SNV using tumor-normal paired samples",
    long_description=open("README.md", encoding="utf-8").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/Xinglab/SomatiX",
    license="SomatiX Research Use License",
    packages=setuptools.find_packages(where="source"),
    package_dir={"": "source"},
    include_package_data=True,
    entry_points={
        "console_scripts": [
            "somatix=somatix.somatix:main",
        ]
    },
    scripts=_scripts,
    classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.12",
        "License :: Other/Proprietary License",
        "Operating System :: POSIX :: Linux",
        "Environment :: Console",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    # Strictly Python 3.12.x
    python_requires=">=3.12,<3.13",
)
