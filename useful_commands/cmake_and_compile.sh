#!/bin/bash

#echo "LD_LIBRARY_PATH=$LD_LIBRARY_PATH"

cmake -DCMAKE_INSTALL_PREFIX=${PWD/%_build}_install ${PWD/%_build} -DGGS_DIR=$GGS_SYS/

make
#make install

rm -fv lib/TestGeometry.gdml
GGSWolowitz -g lib/libTestGeometry.so -gd macros/geo.mac -o lib/TestGeometry.gdml -t gdml

rm -fv lib/TestGeometry_nocalo.gdml
GGSWolowitz -g lib/libTestGeometry_nocalo.so -gd macros/geo.mac -o lib/TestGeometry_nocalo.gdml -t gdml
