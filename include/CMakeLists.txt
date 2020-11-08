set(CMAKE_CXX_FLAGS "${GGS_CXX_FLAGS} ${CMAKE_CXX_FLAGS}")
set_property(DIRECTORY ${CMAKE_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES ${include_dirs})

configure_file(TrCluster.hh TrCluster.hh COPYONLY)
configure_file(TrCluster_linkdef.h TrCluster_linkdef.h COPYONLY)
#configure_file(physics.h physics.h COPYONLY)
#configure_file(absStrip.h absStrip.h COPYONLY)
#configure_file(vector2.h, vector2.h COPYONLY)

install(FILES
        TrCluster.hh
        TrCluster_linkdef.h
        physics.h
        absStrip.h
        vector2.h

        DESTINATION "${INSTALL_INC_DIR}"
        )
