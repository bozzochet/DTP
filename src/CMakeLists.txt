include_directories(${PROJECT_SOURCE_DIR}/include)

####################### Build HerdBaselineGeometry shared library ##################
# Set source

set(TestGeometry_SOURCE
  DetectorConstruction.cc
  DetectorMessenger.cc
)

# Set target
add_library(TestGeometry SHARED ${TestGeometry_SOURCE} )

target_link_libraries(TestGeometry
                      ${Geant4_LIBRARIES}
                      ${GGS_GEO_LIBRARIES}
                     )

set_target_properties(TestGeometry PROPERTIES VERSION 1)

add_library(TestGeometry_nocalo SHARED ${TestGeometry_SOURCE} )
target_compile_definitions(TestGeometry_nocalo PRIVATE _NOCALO_)
target_link_libraries(TestGeometry_nocalo
  ${Geant4_LIBRARIES}
  ${GGS_GEO_LIBRARIES}
)
set_target_properties(TestGeometry_nocalo PROPERTIES VERSION 1)


# Set install rules
install(TARGETS TestGeometry DESTINATION plugins)

find_package(ROOT)

include( ${ROOT_USE_FILE} )
root_generate_dictionary( G__TrCluster TrCluster.hh LINKDEF ${CMAKE_SOURCE_DIR}/include/TrCluster_linkdef.h)

set(SimpleAnalysis_SRC
	SimpleAnalysis1.cpp
	G__TrCluster.cxx
)

add_executable(SimpleAnalysis ${SimpleAnalysis_SRC})
target_link_libraries(SimpleAnalysis
  ${ROOT_LIBRARIES}
  ${GGS_ANALYSIS_LIBRARIES}
)

set(DataAnalysis_SRC
	DataAnalysis.cpp
	G__TrCluster.cxx
)

add_executable(DataAnalysis ${DataAnalysis_SRC})
target_link_libraries(DataAnalysis
  ${ROOT_LIBRARIES}
  ${GGS_ANALYSIS_LIBRARIES}
)

set(DataAnalysisGhosts_SRC
	DataAnalysisGhosts.cpp
	G__TrCluster.cxx
)

add_executable(DataAnalysisGhosts ${DataAnalysisGhosts_SRC})
target_link_libraries(DataAnalysisGhosts
  ${ROOT_LIBRARIES}
  ${GGS_ANALYSIS_LIBRARIES}
)

set(DataAnalysisHough_SRC
	DataAnalysisHough.cpp
	G__TrCluster.cxx
)

add_executable(DataAnalysisHough ${DataAnalysisHough_SRC})
target_link_libraries(DataAnalysisHough
  ${ROOT_LIBRARIES}
  ${GGS_ANALYSIS_LIBRARIES}
)

set(DataAnalysisMulti_SRC
	DataAnalysisMulti.cpp
	G__TrCluster.cxx
)

add_executable(DataAnalysisMulti ${DataAnalysisMulti_SRC})
target_link_libraries(DataAnalysisMulti
  ${ROOT_LIBRARIES}
  ${GGS_ANALYSIS_LIBRARIES}
)



set(alternative_simple_analysis_SRC
        alternative_simple_analysis.cpp
        G__TrCluster.cxx
)

add_executable(alternative_simple_analysis ${alternative_simple_analysis_SRC})
target_link_libraries(alternative_simple_analysis
  ${ROOT_LIBRARIES}
  ${GGS_ANALYSIS_LIBRARIES}
)
