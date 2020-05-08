echo "LD_LIBRARY_PATH=$LD_LIBRARY_PATH"

GGSPenny -g src/libTestGeometry.so -d run.mac -ro test.root
GGSPenny -g src/libTestGeometry_nocalo.so -d run.mac -ro test_nocalo.root
