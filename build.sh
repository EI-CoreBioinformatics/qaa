#!/bin/bash

version=$1

# rm dist/*whl
python setup.py bdist_wheel
# pip install --prefix=/ei/software/testing/qaa/${version}/x86_64 -U dist/qaa-${version}-*.whl
pip install -U dist/qaa-${version}-*.whl
