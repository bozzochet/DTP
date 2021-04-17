#!/bin/bash

#echo "LD_LIBRARY_PATH=$LD_LIBRARY_PATH"

sed 's|/run/beamOn 5000|/run/beamOn 100|' macros/run.mac > macros/run_short.mac

GGSPenny -g lib/libTestGeometry.so -gd macros/geo.mac -d macros/run_short.mac -ro test.root
GGSPenny -g lib/libTestGeometry_nocalo.so -gd macros/geo.mac -d macros/run_short.mac -ro test_nocalo.root
