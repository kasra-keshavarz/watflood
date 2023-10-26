# This file is part of WATFLOOD
# Copyright (C) 1972-2023 Nick Kouwen
#
# This is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3, or (at your option)
# any later version.
#
# This software is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# Under Section 7 of GPL version 3, you are granted additional
# permissions described in the GCC Runtime Library Exception, version
# 3.1, as published by the Free Software Foundation.
#
# You should have received a copy of the GNU General Public License along with
# this program; see the files LICENSE and LICENSE_EXCEPTION respectively.
# If not, see <http://www.gnu.org/licenses/>.
#
# Last modifications:
#    * Kasra Keshavarz <kasra.keshavarz1@ucalgary.ca> 25 Oct 2023
#
# CMAKE build file for WATFLOOD/CHARM
# Options:
# -D CMAKE_BUILD_TYPE=debug
# -D CMAKE_Fortran_COMPILER=gfortran/ifort

# minimum version set to 3.20.0
cmake_minimum_required(VERSION 3.20)

# project name
project(WATFLOOD Fortran)

# Check whether the user has selected a BUILD_TYPE
# two build types available: DEBUG, RELEASE
# default is RELEASE
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
