#!/bin/bash -euo

mkdir -p ${PREFIX}/bin

sed -i.bak 's|setuptools.find_packages(),|setuptools.find_namespace_packages(where="."),|' setup.py
rm -rf *.bak

# Install PAV
${PYTHON} -m pip install --no-deps --no-build-isolation --no-cache-dir . -vvv

# Install Chloe
git clone https://github.com/ian-small/chloe
git clone https://github.com/ian-small/chloe_references

cp -rf chloe  ${PREFIX}/bin
cp -rf chloe_references ${PREFIX}/bin

#julia --project=/Users/chrisjackson/VSC_projects/chloe /Users/chrisjackson/VSC_projects/chloe/chloe.jl annotate --help