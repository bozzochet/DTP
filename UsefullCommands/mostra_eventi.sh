echo "LD_LIBRARY_PATH=$LD_LIBRARY_PATH"

GGSWolowitz -g src/libTestGeometry.so -t vgm -o src/TestGeometry.vgm.root
#GGSLeonard -g src/TestGeometry.gdml -i test.root
GGSLeonard -g src/TestGeometry.vgm.root -i test.root
