#!/usr/bin/env bash

find . -type f -name "._*" -delete
find . -type f -name ".DS_Store" -delete
find . -type f -name "*.py[co]" -delete
find . -type d -name "__pycache__" -delete
find . -type d -name ".pytest_cache" -exec rm -rf {} \;

rm -rf build dist *.spec
