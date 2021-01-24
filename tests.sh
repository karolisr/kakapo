#!/usr/bin/env bash
cd /home/travis/build/karolisr/kakapo/kakapo/utils/c/src
make test
cd /home/travis/build/karolisr/kakapo/tests
pytest
