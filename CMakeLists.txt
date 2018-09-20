cmake_minimum_required(VERSION 3.0)
project(UnigenQA)
set(PROJECT_VERSION 1.0)
set(CMAKE_cpp_STANDARD 11)
#set(CMAKE_BUILD_TYPE debug)

file(GLOB UNIGEN_SOURCES Unigen/*.cxx)
file(GLOB UNIGEN_HEADERS Unigen/*.h)
file(GLOB UNIGEN_HEADERS_DICT RELATIVE ${CMAKE_SOURCE_DIR}/Unigen Unigen/*.h)
file(GLOB QA_SOURCES src/*.cxx)
file(GLOB QA_HEADERS src/*.h)
file(GLOB QA_HEADERS_DICT RELATIVE ${CMAKE_SOURCE_DIR}/src src/*.h)
list(REMOVE_ITEM UNIGEN_HEADERS_DICT UnigenLinkDef.h)

set(HEADERS
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

add_library(UnigenQA SHARED ${QA_SOURCES})
target_link_libraries(UnigenQA ${ROOT_LIBRARIES} Unigen)
add_dependencies(UnigenQA Unigen)

# Create a main program using the library
add_executable(RunUnigenQA src/RunUnigenQA.C)
target_link_libraries(RunUnigenQA UnigenQA ${ROOT_LIBRARIES})

add_executable(RunUnigenQAFilelist src/RunUnigenQAFilelist.C)
target_link_libraries(RunUnigenQAFilelist UnigenQA ${ROOT_LIBRARIES})







include(GenerateExportHeader)


generate_export_header(UnigenQA)
set_property(TARGET UnigenQA PROPERTY VERSION ${PROJECT_VERSION})
set_property(TARGET UnigenQA PROPERTY SOVERSION ${PROJECT_VERSION})
set_property(TARGET UnigenQA PROPERTY INTERFACE_UnigenQA_MAJOR_VERSION 1)
set_property(TARGET UnigenQA APPEND PROPERTY COMPATIBLE_INTERFACE_STRING INTERFACE_UnigenQA_MAJOR_VERSION)

install(TARGETS Unigen UnigenQA
        EXPORT UnigenQATargets
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
        RUNTIME DESTINATION bin
        INCLUDES DESTINATION include
        )

install(
        FILES
        ${HEADERS}
        "${CMAKE_CURRENT_BINARY_DIR}/unigenqa_export.h"
        DESTINATION
        include
        COMPONENT
        Devel
)

#install(
#        FILES
#        "${CMAKE_CURRENT_BINARY_DIR}/${PCM_FILE_NAME}_rdict.pcm"
#        DESTINATION
#        lib
#        OPTIONAL
#)

include(CMakePackageConfigHelpers)
write_basic_package_version_file(
        "${CMAKE_CURRENT_BINARY_DIR}/UnigenQA/UnigenQAConfigVersion.cmake"
        VERSION ${PROJECT_VERSION}
        COMPATIBILITY AnyNewerVersion
)

export(EXPORT UnigenQATargets
        FILE "${CMAKE_CURRENT_BINARY_DIR}/UnigenQA/UnigenQATargets.cmake")

set(ConfigPackageLocation lib/cmake/UnigenQA)

set(UnigenQA_INCLUDE_DIR "include")
set(UnigenQA_LIBRARY_DIR "lib")

configure_package_config_file(
        UnigenQAConfig.cmake.in ${CMAKE_CURRENT_BINARY_DIR}/UnigenQA/UnigenQAConfig.cmake
        INSTALL_DESTINATION ${ConfigPackageLocation}
        PATH_VARS UnigenQA_INCLUDE_DIR UnigenQA_LIBRARY_DIR
        #  [NO_SET_AND_CHECK_MACRO]
        #  [NO_CHECK_REQUIRED_COMPONENTS_MACRO]
        #  [INSTALL_PREFIX <path>]
)

install(EXPORT UnigenQATargets
        FILE
        UnigenQATargets.cmake
        DESTINATION
        ${ConfigPackageLocation}
        )

install(FILES
        "${CMAKE_CURRENT_BINARY_DIR}/UnigenQA/UnigenQAConfig.cmake"
        "${CMAKE_CURRENT_BINARY_DIR}/UnigenQA/UnigenQAConfigVersion.cmake"
        DESTINATION
        ${ConfigPackageLocation}
        COMPONENT
        Devel
        )



