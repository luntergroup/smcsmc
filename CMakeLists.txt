cmake_minimum_required(VERSION 3.1)

set(Boost_DEBUG 1) 
set (CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake)

project(smcsmc)

set(VERSION 0.1-beta)

set(CMAKE_CXX_STANDARD 11)

set(SCRM_SRC src/scrm/src/node.cc 
             src/scrm/src/node.h
             src/scrm/src/node_container.cc
             src/scrm/src/node_container.h
             src/scrm/src/time_interval.cc
             src/scrm/src/time_interval.h
             src/scrm/src/tree_point.cc
             src/scrm/src/tree_point.h
             src/scrm/src/forest.cc
             src/scrm/src/forest-debug.cc
             src/scrm/src/random/mersenne_twister.h
             src/scrm/src/random/random_generator.cc
             src/scrm/src/random/mersenne_twister.cc
             src/scrm/src/random/random_generator.h
             src/scrm/src/random/fastfunc.cc
             src/scrm/src/random/fastfunc.h
             src/scrm/src/param.cc
             src/scrm/src/param.h
             src/scrm/src/model.cc
             src/scrm/src/model.h
             src/scrm/src/summary_statistics/newick_tree.cc
             src/scrm/src/summary_statistics/frequency_spectrum.cc
             src/scrm/src/summary_statistics/seg_sites.cc
             src/scrm/src/summary_statistics/tmrca.cc
             src/scrm/src/summary_statistics/oriented_forest.cc)

set(SMCSMC_SRC src/coalevent.cpp
               src/particle.cpp
               src/particleContainer.cpp
               src/count.cpp
               src/count-semi-debug.cpp
               src/pattern.cpp
               src/pfparam.cpp
               src/segdata.cpp
               src/arena.cpp
               src/model_summary.cpp)

execute_process(COMMAND git submodule update --init --recursive)

execute_process(COMMAND bash -c "git submodule status | grep scrm | cut -d' ' -f4 | sed 's/[)(]//g' | tr -d '\n'" OUTPUT_VARIABLE SCRMVERSION)
execute_process(COMMAND bash -c "git describe --tags | tr -d '\n'" OUTPUT_VARIABLE SMCSMCVERSION)
execute_process(COMMAND bash -c "date -u | sed -e 's/ /-/g' | tr -d '\n'" OUTPUT_VARIABLE COMPILEDATE)

string(CONCAT CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}"
			      " -DNDEBUG "
                              "-Wall -Wno-unknown-pragmas -pedantic -g -O3 -DLANE_MAX_SIZE=50 "
                              "-Wno-deprecated-declarations "                     # avoid complaints about auto_ptr
                              "-DSMCSMCVERSION=\\\"${SMCSMCVERSION}\\\" "
                              "-DSCRMVERSION=\\\"${SCRMVERSION}\\\" "
                              "-DCOMPILEDATE=\\\"${COMPILEDATE}\\\" "
                              "-DVERSION=\\\"${VERSION}\\\" ")

find_package(Boost 1.67.0 COMPONENTS iostreams regex)
find_package(ZLIB)
include_directories(${Boost_INCLUDE_DIRS})
include_directories(${CONDA_PREFIX}/include)

add_executable(smcsmc ${SMCSMC_SRC} ${SCRM_SRC} src/smcsmc.cpp)
add_executable(scrm ${SCRM_SRC} src/scrm/src/scrm.cc)

target_link_libraries(smcsmc ${Boost_LIBRARIES})
message(STATUS "Boost_LIBRARIES: ${Boost_LIBRARIES}")
target_link_libraries(smcsmc ${ZLIB_LIBRARIES})
target_include_directories(smcsmc PRIVATE src/scrm/src src/scrm/src/random)
#SET(CMAKE_INSTALL_RPATH "${CMAKE_INSaTALL_PREFIX}/lib")

install(TARGETS smcsmc scrm
	DESTINATION bin)
