echo "LD_LIBRARY_PATH=$LD_LIBRARY_PATH"

GGSPenny -g src/libTestGeometry.dylib -d run.mac -ro test.root
GGSPenny -g src/libTestGeometry_nocalo.dylib -d run.mac -ro test_nocalo.root
