echo "LD_LIBRARY_PATH=$LD_LIBRARY_PATH"

cmake ../myggsrepo -DGGS_DIR=../../GGSSoftware_install/
make
rm -v src/TestGeometry.gdml
GGSWolowitz -g src/libTestGeometry.dylib -o src/TestGeometry.gdml -t gdml
GGSLeonard -G src/TestGeometry.gdml
