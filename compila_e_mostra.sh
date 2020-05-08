echo "LD_LIBRARY_PATH=$LD_LIBRARY_PATH"

cmake ${PWD/%_build} -DGGS_DIR=$GGS_SYS/
make
rm -v src/TestGeometry.gdml
GGSWolowitz -g src/libTestGeometry.so -o src/TestGeometry.gdml -t gdml
#GGSLeonard -G src/TestGeometry.gdml
rm -v src/TestGeometry_nocalo.gdml
GGSWolowitz -g src/libTestGeometry_nocalo.so -o src/TestGeometry_nocalo.gdml -t gdml
#GGSLeonard -G src/TestGeometry_nocalo.gdml
