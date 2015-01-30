cmake_minimum_required (VERSION 2.8)

set(CMAKE_C_COMPILER "icc")
set(CMAKE_CXX_COMPILER "icpc")
list( APPEND CMAKE_CXX_FLAGS "-std=c++0x ${CMAKE_CXX_FLAGS}")
SET(CMAKE_INSTALL_PREFIX /Users/brbass/Documents/college3/research/spn_transport/manatee/bin)

project(Manatee)

enable_testing()

find_package(MPI REQUIRED)
find_package(Trilinos PATHS ${TRILINOS_PATH}/lib/cmake/Trilinos ${TRILINOS_PATH} REQUIRED)

IF(Trilinos_FOUND)
        MESSAGE("\nFound Trilinos!  Here are the details: ")
        MESSAGE("   Trilinos_DIR = ${Trilinos_DIR}")
        MESSAGE("   Trilinos_VERSION = ${Trilinos_VERSION}")
        MESSAGE("   Trilinos_PACKAGE_LIST = ${Trilinos_PACKAGE_LIST}")
        MESSAGE("   Trilinos_LIBRARIES = ${Trilinos_LIBRARIES}")
        MESSAGE("   Trilinos_INCLUDE_DIRS = ${Trilinos_INCLUDE_DIRS}")
        MESSAGE("   Trilinos_TPL_LIST = ${Trilinos_TPL_LIST}")
        MESSAGE("   Trilinos_TPL_INCLUDE_DIRS = ${Trilinos_TPL_INCLUDE_DIRS}")
        MESSAGE("   Trilinos_TPL_LIBRARIES = ${Trilinos_TPL_LIBRARIES}")
        MESSAGE("   Trilinos_BUILD_SHARED_LIBS = ${Trilinos_BUILD_SHARED_LIBS}")
        MESSAGE("   Trilinos_CXX_COMPILER = ${Trilinos_CXX_COMPILER}")
        MESSAGE("   Trilinos_C_COMPILER = ${Trilinos_C_COMPILER}")
        MESSAGE("   Trilinos_Fortran_COMPILER = ${Trilinos_Fortran_COMPILER}")
        MESSAGE("   Trilinos_CXX_COMPILER_FLAGS = ${Trilinos_CXX_COMPILER_FLAGS}")
        MESSAGE("   Trilinos_C_COMPILER_FLAGS = ${Trilinos_C_COMPILER_FLAGS}")
        MESSAGE("   Trilinos_Fortran_COMPILER_FLAGS = ${Trilinos_Fortran_COMPILER_FLAGS}")
        MESSAGE("   Trilinos_LINKER = ${Trilinos_LINKER}")
        MESSAGE("   Trilinos_EXTRA_LD_FLAGS = ${Trilinos_EXTRA_LD_FLAGS}")
        MESSAGE("   Trilinos_AR = ${Trilinos_AR}")
        MESSAGE("End of Trilinos details\n")
ELSE()
        MESSAGE(FATAL_ERROR "Could not find Trilinos!")
ENDIF()

set(manatee_src src/data/Data.cc src/mesh/Finite_Element.cc src/mesh/Mesh.cc src/neutronics/Neutronics.cc src/parser/Parser.cc src/transport_model/SPn_Transport.cc)
set(manatee_main src/manatee/Manatee.cc)
set(manatee_test test/tst_Manatee.cc)

include_directories(${MPI_INCLUDE_PATH})
include_directories(${Trilinos_INCLUDE_DIRS} ${Trilinos_TPL_INCLUDE_DIRS})
include_directories(src/data src/mesh src/neutronics src/parser src/transport_model)

link_directories(${Trilinos_LIBRARY_DIRS} ${Trilinos_TPL_LIBRARY_DIRS})

add_executable(manatee ${manatee_src} ${manatee_main})
add_executable(tst_manatee ${manatee_src} ${manatee_test})

add_test(test_manateee tst_manatee)

set(CMAKE_VERBOSE_MAKEFILE on)
set(CMAKE_BUILD_TYPE Release)

target_link_libraries(manatee ${MPI_LIBRARIES})
target_link_libraries(tst_manatee ${MPI_LIBRARIES})
target_link_libraries(manatee ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES})
target_link_libraries(tst_manatee ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES})

install(TARGETS manatee DESTINATION .)
