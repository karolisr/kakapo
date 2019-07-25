#!/bin/bash

find . -type f -name ".DS_Store" -delete
find . -type f -name "*.py[co]" -delete
find . -type d -name "__pycache__" -delete

rm -rf build dist *.spec
