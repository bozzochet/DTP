set(CMAKE_CXX_FLAGS "${GGS_CXX_FLAGS} ${CMAKE_CXX_FLAGS}")
set_property(DIRECTORY ${CMAKE_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES ${include_dirs})

configure_file(TrCluster.hh TrCluster.hh COPYONLY)
configure_file(TrCluster_linkdef.h TrCluster_linkdef.h COPYONLY)
configure_file(Geometry.h Geometry.h COPYONLY)
configure_file(measure.h measure.h COPYONLY)
configure_file(physics.h physics.h COPYONLY)
configure_file(progress.h progress.h COPYONLY)
configure_file(TimeSegm.h TimeSegm.h COPYONLY)
configure_file(TimeSim.h TimeSim.h COPYONLY)
configure_file(vector2.h vector2.h COPYONLY)


install(DIRECTORY . DESTINATION include
       USE_SOURCE_PERMISSIONS
       PATTERN *~ EXCLUDE
       PATTERN "DEBUG.h" EXCLUDE
       PATTERN "CMakeLists.txt" EXCLUDE
       )
