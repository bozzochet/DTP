#!/bin/bash

rm -fv lib/TestGeometry.vgm.root
GGSWolowitz -g lib/libTestGeometry.so -t vgm -o lib/TestGeometry.vgm.root
#GGSLeonard -g lib/TestGeometry.gdml -i test.root
GGSLeonard -g lib/TestGeometry.vgm.root -i test.root

rm -fv lib/TestGeometry_nocalo.vgm.root
GGSWolowitz -g lib/libTestGeometry_nocalo.so -t vgm -o lib/TestGeometry_nocalo.vgm.root
#GGSLeonard -g lib/TestGeometry_nocalo.gdml -i test_nocalo.root
GGSLeonard -g lib/TestGeometry_nocalo.vgm.root -i test_nocalo.root
