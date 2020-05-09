# DTP

## Simulation with GGS

GGS Simulation of a simple detector with the silicon detectors with timing capabilities

Ingredients:

- macros/rin.mac: datacard, a la Geant, to set the simulation parameters
<!-- - macros/geo.mac: datacard for the parametric geometry -->
- {src,include}/DetectorConstruction.{cc,hh}: definition, a la Geant, of the geometry
<!-- - Analysis/Analysis.C: ROOT macro to read the GGS output file -->

Typical commands:

##### - compilation and installation:
```
cd <build_path> 
cmake -DCMAKE_INSTALL_PREFIX=<installation_path> <source_path> -DGGS_DIR=$GGS_SYS/
make
make install
```

The `make install` is not required and even not suggested. All the pieces are both copied in the build directory and installed: everything should work even in the build directory.

A working example can be found in `useful_commands/cmake_and_compile.sh` that includes also the creation of the `gdml` version of the geometry:
```
GGSWolowitz -g lib/libTestGeometry.so -o lib/TestGeometry.gdml -t gdml
```

##### - simulation:
```
GGSPenny -g lib/libTestGeometry.so -d macros/run.mac -ro test.root
```

A working example can be found in `useful_commands/run.sh`

- visualization of the geometry and of the simulated events:
```
GGSLeonard -g lib/TestGeometry.gdml
```
```
GSLeonard -g lib/TestGeometry.gdml -i test.root
```

Working examples can be found in `useful_commands/show_geometry.sh` and  `useful_commands/show_events.sh`

- conversion from parametric geometry to GDML
```
GGSWolowitz -g lib/libTestGeometry.so -o lib/TestGeometry.gdml -t gdml
```

##### - conversion from parametric geometry to VGM (http://ivana.home.cern.ch/ivana/VGM.html)
```
GGSWolowitz -g lib/libTestGeometry.so -t vgm -o lib/TestGeometry.vgm.root
```

## Analysis (to be written correctly)

- conversion from GGS output to plain ROOT file:
```
root [0] .L Analysis/Analysis.C
root [1] SimpleAnalysis("GGSRootOutput.root", "anaOut.root")
