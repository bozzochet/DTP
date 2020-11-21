set(CMAKE_CXX_FLAGS "${GGS_CXX_FLAGS} ${CMAKE_CXX_FLAGS}")
set_property(DIRECTORY ${CMAKE_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES ${include_dirs})

configure_file(run.mac run.mac COPYONLY)

configure_file(e_10.mac e_10.mac COPYONLY)
configure_file(e_100.mac e_100.mac COPYONLY)
configure_file(e_1000.mac e_1000.mac COPYONLY)
configure_file(e_10000.mac e_10000.mac COPYONLY)

configure_file(p_10.mac p_10.mac COPYONLY)
configure_file(p_100.mac p_100.mac COPYONLY)
configure_file(p_1000.mac p_1000.mac COPYONLY)
configure_file(p_10000.mac p_10000.mac COPYONLY)

configure_file(geo.mac geo.mac COPYONLY)
configure_file(ToySegmentation.C ToySegmentation.C COPYONLY)
configure_file(ToyTiming.C ToyTiming.C COPYONLY)
configure_file(ElectronProton.C ElectronProton.C COPYONLY)

install(DIRECTORY . DESTINATION macros
       USE_SOURCE_PERMISSIONS
       PATTERN *~ EXCLUDE
       PATTERN "CMakeLists.txt" EXCLUDE
       )
