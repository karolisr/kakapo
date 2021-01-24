#!/usr/bin/env bash
cd /home/travis/build/karolisr/kakapo/kakapo/utils/c/src
make deep-clean
make test
make
cd /home/travis/build/karolisr/kakapo/tests
pytest
