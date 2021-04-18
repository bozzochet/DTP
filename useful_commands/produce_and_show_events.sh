#!/bin/bash

#echo "LD_LIBRARY_PATH=$LD_LIBRARY_PATH"

sed 's|/run/beamOn 5000||' macros/run.mac > macros/run_short.mac

GGSPenny -g lib/libTestGeometry.so -gd macros/geo.mac -d macros/run_short.mac -X
