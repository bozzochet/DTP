echo "LD_LIBRARY_PATH=$LD_LIBRARY_PATH"

cmake ../myggsrepo -DGGS_DIR=../../GGSSoftware_install/
make
rm -v src/TestGeometry.gdml
GGSWolowitz -g src/libTestGeometry.dylib -o src/TestGeometry.gdml -t gdml
#GGSLeonard -G src/TestGeometry.gdml
rm -v src/TestGeometry_nocalo.gdml
GGSWolowitz -g src/libTestGeometry_nocalo.dylib -o src/TestGeometry_nocalo.gdml -t gdml
#GGSLeonard -G src/TestGeometry_nocalo.gdml
