#!/bin/bash

#echo "LD_LIBRARY_PATH=$LD_LIBRARY_PATH"

GGSPenny -g lib/libTestGeometry.so -d macros/run.mac -ro test.root
GGSPenny -g lib/libTestGeometry_nocalo.so -d macros/run.mac -ro test_nocalo.root
