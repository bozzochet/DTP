project(TestGeometry)

find_package(GGS 2.3.1 REQUIRED)
set(CMAKE_CXX_FLAGS "${GGS_CXX_FLAGS} ${CMAKE_CXX_FLAGS}")
# Remove /usr/include from list of include folders.
# This folder is added by Geant4 through GGS, and it causes this error:
#   #error "You need a ISO C conforming compiler to use the glibc headers"
# when generating Root dictionaries in a system with a recent version of glibc
# (tested with glibc 2.23).
get_property(include_dirs DIRECTORY ${CMAKE_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES)
list(REMOVE_ITEM include_dirs "/usr/include")
set_property(DIRECTORY ${CMAKE_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES ${include_dirs})


add_subdirectory(src)
add_subdirectory(ManualMakefile)
configure_file(compila_e_mostra.sh compila_e_mostra.sh COPYONLY)
configure_file(runna.sh runna.sh COPYONLY)
configure_file(run.mac run.mac COPYONLY)
