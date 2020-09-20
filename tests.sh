#!/usr/bin/env bash
cd kakapo/utils/c
make deep-clean
make
cd ../../../
cd tests
pytest
