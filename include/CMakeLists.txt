set(CMAKE_CXX_FLAGS "${GGS_CXX_FLAGS} ${CMAKE_CXX_FLAGS}")
set_property(DIRECTORY ${CMAKE_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES ${include_dirs})

configure_file(TrCluster.hh TrCluster.hh COPYONLY)
configure_file(TrCluster_linkdef.h TrCluster_linkdef.h COPYONLY)

install(FILES TrCluster.hh DESTINATION "${INSTALL_INC_DIR}")
install(FILES TrCluster_linkdef.h DESTINATION "${INSTALL_INC_DIR}")
