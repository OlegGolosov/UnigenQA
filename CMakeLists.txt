cmake_minimum_required(VERSION 2.4)
project(UnigenQA)
set(PROJECT_VERSION 1.0)

file(GLOB UNIGEN_SOURCES Unigen/*.cxx)
file(GLOB UNIGEN_HEADERS Unigen/*.h)
file(GLOB UNIGEN_HEADERS_DICT RELATIVE ${CMAKE_SOURCE_DIR}/Unigen Unigen/*.h)
list(REMOVE_ITEM UNIGEN_HEADERS_DICT UnigenLinkDef.h)

file(GLOB QA_SOURCES src/*.cxx)
file(GLOB QA_HEADERS src/*.h)
file(GLOB QA_HEADERS_DICT RELATIVE ${CMAKE_SOURCE_DIR}/src src/*.h)
list(REMOVE_ITEM QA_HEADERS_DICT UnigenQALinkDef.h)

set (HEADERS 
	${UNIGEN_HEADERS}
	${QA_HEADERS}
)

set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")

# Telling CMake where to find the ROOT installation.
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
list(APPEND CMAKE_PREFIX_PATH "${CMAKE_SOURCE_DIR}/lib")

#---Locate the ROOT package and defines a number of variables (e.g. ROOT_INCLUDE_DIRS)
find_package(ROOT REQUIRED COMPONENTS Core MathCore RIO)
#find_package(ROOT REQUIRED COMPONENTS Core MathCore MathMore RIO Hist Tree Net TreePlayer)

# Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
include(${ROOT_USE_FILE})
message("${ROOT_USE_FILE}")
add_definitions(${ROOT_CXX_FLAGS})

include_directories(${CMAKE_SOURCE_DIR} ${ROOT_INCLUDE_DIRS} Unigen src)

ROOT_GENERATE_DICTIONARY(G__Unigen ${UNIGEN_HEADERS_DICT} LINKDEF Unigen/UnigenLinkDef.h)
add_library(Unigen SHARED ${UNIGEN_SOURCES} G__Unigen)
target_link_libraries(Unigen ${ROOT_LIBRARIES})

ROOT_GENERATE_DICTIONARY(G__UnigenQA ${QA_HEADERS_DICT} LINKDEF src/UnigenQALinkDef.h)
add_library(UnigenQA SHARED ${QA_SOURCES})
target_link_libraries(UnigenQA Unigen ${ROOT_LIBRARIES})
add_dependencies(UnigenQA Unigen)

# Create a main program using the library
add_executable(RunUnigenQA macro/RunUnigenQA.C)
target_link_libraries(RunUnigenQA UnigenQA Unigen ${ROOT_LIBRARIES})
target_include_directories(RunUnigenQA PUBLIC)
