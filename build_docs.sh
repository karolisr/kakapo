#!/usr/bin/env bash

sphinx-apidoc -A "Karolis Ramanauskas" -H "Kakapo" -V "1.0.0" -F -o docs \
./ \
tests/ \
tests.py \
setup.py \
kakapo/__main__.py

##################################################################################
# This is for BSD variant of sed, which is found on Mac OS X
# the two double quotes at the beginning of the comment are not required otherwise
sed -i "" "s|# import os|import os|g" docs/conf.py
sed -i "" "s|# import sys|import sys|g" docs/conf.py
sed -i "" "s|# sys.path.insert|sys.path.insert|g" docs/conf.py
##################################################################################

# sed -i "s|# import os|import os|g" docs/conf.py
# sed -i "s|# import sys|import sys|g" docs/conf.py
# sed -i "s|# sys.path.insert|sys.path.insert|g" docs/conf.py

cd docs
make html
# make latex
# cd _build/latex
# make
