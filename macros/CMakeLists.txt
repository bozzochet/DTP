set(CMAKE_CXX_FLAGS "${GGS_CXX_FLAGS} ${CMAKE_CXX_FLAGS}")
set_property(DIRECTORY ${CMAKE_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES ${include_dirs})

configure_file(run.mac run.mac COPYONLY)
configure_file(ToySegmentation.C ToySegmentation.C COPYONLY)

install(DIRECTORY . DESTINATION macros
       USE_SOURCE_PERMISSIONS
       PATTERN *~ EXCLUDE
       PATTERN "CMakeLists.txt" EXCLUDE
       )
