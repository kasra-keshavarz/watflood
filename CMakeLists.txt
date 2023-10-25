cmake_minimum_required(VERSION 3.16)

project(WATFLOOD Fortran)

# Check whether the user has selected a BUILD_TYPE
IF ("${CMAKE_BUILD_TYPE}" STREQUAL "DEBUG")
    set(CMAKE_BUILD_TYPE "DEBUG")
ELSE()
    set(CMAKE_BUILD_TYPE "RELEASE")
ENDIF()
message(STATUS "Build type is set to ${CMAKE_BUILD_TYPE}")
message(STATUS "Change build type with -DCMAKE_BUILD_TYPE=DEBUG/RELEASE")

set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/include)

include_directories(${CMAKE_Fortran_MODULE_DIRECTORY})
file(MAKE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY})

set(CMAKE_INSTALL_PREFIX ${PROJECT_BINARY_DIR})

# if DEBUG requested, set the various debug options.
IF(${CMAKE_BUILD_TYPE} MATCHES "DEBUG")
    set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -O0")
    IF(${CMAKE_Fortran_COMPILER_ID} MATCHES "GNU")
        set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -Wall -fcheck=all -fbacktrace")
    ELSEIF(${CMAKE_Fortran_COMPILER_ID} MATCHES "Intel")
        set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -warn all -check all -traceback")
    ELSE()
        message(WARNING "Unable to determine Compiler ID: ${CMAKE_Fortran_COMPILER_ID}")
    ENDIF(${CMAKE_Fortran_COMPILER_ID} MATCHES "GNU")
    message(STATUS "Debug flags for ${CMAKE_Fortran_COMPILER_ID}: ${CMAKE_Fortran_FLAGS_DEBUG}")
ENDIF(${CMAKE_BUILD_TYPE} MATCHES "DEBUG")

add_subdirectory(src)
