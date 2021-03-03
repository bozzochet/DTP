


############## NOT WORKING ##########################



set(CMAKE_CXX_FLAGS "${GGS_CXX_FLAGS} ${CMAKE_CXX_FLAGS}")
include_directories(${PROJECT_SOURCE_DIR}/include)

#################

set(weightfield_SRC
	weightfield2_signal.cpp
)

add_executable(weightfield ${weightfield_SRC})
target_link_libraries(weightfield
  ${ROOT_LIBRARIES}
)

set_target_properties(weightfield PROPERTIES
				     RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/misc"
				     )

install(TARGETS DataAnalysis DESTINATION "misc")

#################

set(SimpleTimeAnalysis_SRC
	SimpleTimeAnalysis.cpp
)

add_executable(SimpleTimeAnalysis ${SimpleTimeAnalysis_SRC})
target_link_libraries(SimpleTimeAnalysis
  ${ROOT_LIBRARIES}
)

set_target_properties(SimpleTimeAnalysis PROPERTIES
				     RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/misc"
				     )

install(TARGETS SimpleTimeAnalysis DESTINATION "misc")
