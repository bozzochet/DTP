#!/bin/bash

#echo "LD_LIBRARY_PATH=$LD_LIBRARY_PATH"

GGSPenny -g lib/libTestGeometry.so -gd macros/geo.mac -d macros/run.mac -ro test.root
GGSPenny -g lib/libTestGeometry_nocalo.so -gd macros/geo.mac -d macros/run.mac -ro test_nocalo.root
