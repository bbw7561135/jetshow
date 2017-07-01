#!/bin/sh

cd cmake-build-debug
cmake ..
make jetshow_test
./jetshow_test