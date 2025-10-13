#!/bin/bash
git init

BRANCH=main
echo "Fetching $BRANCH"
git remote add -f origin git@github.com:indra-ehep/hgcal-tpg-fe.git
git config core.sparseCheckout true
echo "TPGStage2Emulation/*.hh" > .git/info/sparse-checkout
echo "inc/*.hh" > .git/info/sparse-checkout
git checkout main
mv TPGStage2Emulation/* .
mv inc/* .
