#!/bin/sh

git stash -q --keep-index
#bash ./run_tests.sh
RESULT=$?
git stash pop -q
[ $RESULT -ne 0 ] && exit 1
exit 0
