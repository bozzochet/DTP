include_directories(${PROJECT_SOURCE_DIR}/include)

####################### Build Geometry shared library ##################
# Set source

set(TestGeometry_SOURCE
  DetectorConstruction.cc
)

# Set target TestGeometry and TestGeometry_nocalo
add_library(TestGeometry SHARED ${TestGeometry_SOURCE} )
add_library(TestGeometry_nocalo SHARED ${TestGeometry_SOURCE} )

target_link_libraries(TestGeometry
		      ${Geant4_LIBRARIES}
		      ${GGS_GEO_LIBRARIES}
		      )

target_compile_definitions(TestGeometry_nocalo PRIVATE _NOCALO_)
target_link_libraries(TestGeometry_nocalo
		      ${Geant4_LIBRARIES}
		      ${GGS_GEO_LIBRARIES}
		      )

set_target_properties(TestGeometry PROPERTIES
				   VERSION 1
				   LIBRARY_OUTPUT_DIRECTORY "${LIB_DIR}"
				   )
set_target_properties(TestGeometry_nocalo PROPERTIES
					  VERSION 1
					  LIBRARY_OUTPUT_DIRECTORY "${LIB_DIR}"
					  )
if (APPLE)
#    set_property(TARGET TestGeometry PROPERTY PREFIX "lib")
#    set_property(TARGET TestGeometry PROPERTY OUTPUT_NAME "TestGeometry.so")
    set_property(TARGET TestGeometry PROPERTY SUFFIX ".so")
    set_property(TARGET TestGeometry_nocalo PROPERTY SUFFIX ".so")
#    set_property(TARGET TestGeometry PROPERTY SOVERSION "32.1.2.0")
endif()

# Set install rules
install(TARGETS TestGeometry DESTINATION "${INSTALL_LIB_DIR}")
install(TARGETS TestGeometry_nocalo DESTINATION "${INSTALL_LIB_DIR}")

#############################################################################
######################### Other stuff #######################################

find_package(ROOT CONFIG REQUIRED)

include( ${ROOT_USE_FILE} )
include("${ROOT_DIR}/modules/RootNewMacros.cmake")
root_generate_dictionary(G__TrCluster
			TrCluster.hh
			LINKDEF ${CMAKE_SOURCE_DIR}/include/TrCluster_linkdef.h
			MODULE TrCluster
			)

#Install the ROOT 6 PCM and rootmap files.
if (${ROOT_VERSION} VERSION_GREATER "6.0")
   install(
      FILES
      ${CMAKE_CURRENT_BINARY_DIR}/libTrCluster_rdict.pcm
      ${CMAKE_CURRENT_BINARY_DIR}/libTrCluster.rootmap
      DESTINATION "${INSTALL_LIB_DIR}")
endif (${ROOT_VERSION} VERSION_GREATER "6.0")

######################### EXTRA LIBS #########################

set(TIME_SIM_SRC TimeSim.cpp)
add_library(TimeSim SHARED ${TIME_SIM_SRC})
set_target_properties(TimeSim PROPERTIES
  LIBRARY_OUTPUT_DIRECTORY "${LIB_DIR}"
)
install(TARGETS TimeSim DESTINATION "${INSTALL_LIB_DIR}")

#################

set(POS_SIM_SRC PosSim.cpp)
add_library(PosSim SHARED ${POS_SIM_SRC})
set_target_properties(PosSim PROPERTIES
  LIBRARY_OUTPUT_DIRECTORY "${LIB_DIR}"
)
install(TARGETS PosSim DESTINATION "${INSTALL_LIB_DIR}")

#################

set(TIME_SEGM_SRC TimeSegm.cpp)
add_library(TimeSegm SHARED ${TIME_SEGM_SRC})
set_target_properties(TimeSegm PROPERTIES
  LIBRARY_OUTPUT_DIRECTORY "${LIB_DIR}"
)
install(TARGETS TimeSegm DESTINATION "${INSTALL_LIB_DIR}")

#################

set(TIME_LIBRARIES TimeSim TimeSegm)
set(POS_LIBRARIES PosSim)

#################

set(Digitization_SRC
	Digitization.cpp
  progress.cpp
	G__TrCluster.cxx
)

add_executable(Digitization ${Digitization_SRC})
target_link_libraries(Digitization
  ${ROOT_LIBRARIES}
  ${GGS_ANALYSIS_LIBRARIES}
  ${TIME_LIBRARIES}
  ${POS_LIBRARIES}
)

set_target_properties(Digitization PROPERTIES
				     RUNTIME_OUTPUT_DIRECTORY "${BIN_DIR}"
				     )

install(TARGETS Digitization DESTINATION "${INSTALL_BIN_DIR}")

#################

set(DataAnalysis_SRC
	DataAnalysis.cpp
  progress.cpp
	G__TrCluster.cxx
)

add_executable(DataAnalysis ${DataAnalysis_SRC})
target_link_libraries(DataAnalysis
  ${ROOT_LIBRARIES}
  ${GGS_ANALYSIS_LIBRARIES}
)

set_target_properties(DataAnalysis PROPERTIES
				     RUNTIME_OUTPUT_DIRECTORY "${BIN_DIR}"
				     )

install(TARGETS DataAnalysis DESTINATION "${INSTALL_BIN_DIR}")

#################
