#   Created by Denis-Michael Lux on 05. November 2015.
#
#   This file is part of PDSOFTlib.
#
#   PDSOFTlib is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   PDSOFTlib is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with PDSOFT.  If not, see <http://www.gnu.org/licenses/>.

CMAKE_MINIMUM_REQUIRED(VERSION 2.8 FATAL_ERROR)

SET(PDSOFT_MAJOR 1)
SET(PDSOFT_MINOR 0)
SET(PDSOFT_PATCH 0)

SET(PDSOFT_VERSION ${PDSOFT_MAJOR}.${PDSOFT_MINOR}.${PDSOFT_PATCH})
MESSAGE(STATUS "Configuring PDSOFT ${PDSOFT_VERSION}")

SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake_aux/")

PROJECT(PDSOFT CXX)

MESSAGE(STATUS "")
MESSAGE(STATUS "CMAKE_SYSTEM_NAME          = ${CMAKE_SYSTEM_NAME}"         )
MESSAGE(STATUS "CMAKE_CXX_COMPILER_ID      = ${CMAKE_CXX_COMPILER_ID}"     )
MESSAGE(STATUS "CMAKE_CXX_COMPILER_VERSION = ${CMAKE_CXX_COMPILER_VERSION}")
MESSAGE(STATUS "CMAKE_COMPILER_IS_GNUCXX   = ${CMAKE_COMPILER_IS_GNUCXX}"  )
MESSAGE(STATUS "")

MESSAGE(STATUS "*** Try to find FFTW")

# Including FFTW
INCLUDE(PDSOFTlib_FindFFTW)
MESSAGE(STATUS "~> FFTW_FOUND = ${FFTW_FOUND}")

IF(FFTW_FOUND)
    SET(PDSOFT_LIBS ${PDSOFT_LIBS} ${FFTW_LIB})
    SET(PDSOFT_INCLUDE_DIRS ${PDSOFT_INCLUDE_DIRS} ${FFTW_INCLUDE_DIR})
ENDIF()

MESSAGE(STATUS "")
MESSAGE(STATUS "*** PDSOFT wrapper library will use the following libraries:")
MESSAGE(STATUS "*** PDSOFT_LIBS          = ${PDSOFT_LIBS}"                  )
MESSAGE(STATUS "*** PDSOFT_INCLUDE_DIRS  = ${PDSOFT_INCLUDE_DIRS}"          )
MESSAGE(STATUS "")

IF(CMAKE_COMPILER_IS_GNUCXX)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
    MESSAGE(STATUS "Detected gcc. Added '-O3' to compiler flags")
ENDIF()

IF(DEFINED CMAKE_CXX_COMPILER_ID AND DEFINED CMAKE_CXX_COMPILER_VERSION)
    IF(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND NOT ${CMAKE_CXX_COMPILER_VERSION} VERSION_LESS 4.8.1)
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")

        MESSAGE(STATUS "Detected gcc 4.8.1 or later. Added '-std=c++11 -fopenmp' to compiler flags")
        MESSAGE(STATUS "OpenMP support for PDSOFT enabled.")
    ENDIF()
ENDIF()

INCLUDE_DIRECTORIES( ${CMAKE_CURRENT_SOURCE_DIR}/include                )
INCLUDE_DIRECTORIES( ${CMAKE_CURRENT_SOURCE_DIR}/include/PDSOFT_headers )
INCLUDE_DIRECTORIES( ${PDSOFT_INCLUDE_DIRS}                             )

IF(APPLE)
    IF("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
        MESSAGE(STATUS "Clang compiler on MacOS X detected. Added '-stdlib=libc++ -std=c++11 -O3' to compiler flags.")
    ENDIF()
ENDIF()

MESSAGE(STATUS "")
MESSAGE(STATUS "CMAKE_CXX_FLAGS           = ${CMAKE_CXX_FLAGS}"          )
MESSAGE(STATUS "CMAKE_SHARED_LINKER_FLAGS = ${CMAKE_SHARED_LINKER_FLAGS}")
MESSAGE(STATUS "CMAKE_REQUIRED_INCLUDES   = ${CMAKE_REQUIRED_INCLUDES}"  )
MESSAGE(STATUS "")

SET(PDSOFT_DEBUG           1 CACHE BOOL "Show debugging information in console if they occure at execution time.")
SET(PDSOFT_SHOW_WARNINGS   1 CACHE BOOL "Show warning messages in console if they occure at execution time.")
SET(PDSOFT_SHOW_ERRORS     1 CACHE BOOL "show error messages in console if they occure at execution time.")

MESSAGE(STATUS "")
MESSAGE(STATUS "*** Configure the compiler_config.hpp.in file. Filling")
MESSAGE(STATUS "*** CMake variable contents into the compiler_config.hpp")
MESSAGE(STATUS "*** header file.")
MESSAGE(STATUS "")

CONFIGURE_FILE(
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake_aux/compiler_config.hpp.in
    ${CMAKE_CURRENT_SOURCE_DIR}/include/PDSOFTlib_headers/compiler_config.hpp
)

FILE(GLOB SRCS ${SRCS} ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp                       )
FILE(GLOB SRCS ${SRCS} ${CMAKE_CURRENT_SOURCE_DIR}/include/PDSOFTlib_headers/*.hpp )
FILE(GLOB SRCS ${SRCS} ${CMAKE_CURRENT_SOURCE_DIR}/include/PDSOFTlib               )

ADD_LIBRARY( PDSOFT SHARED ${SRCS} )
TARGET_LINK_LIBRARIES( PDSOFT ${PDSOFT_LIBS} )
SET_TARGET_PROPERTIES( PDSOFT PROPERTIES VERSION ${PDSOFT_VERSION} SOVERSION ${PDSOFT_MAJOR} )

ADD_SUBDIRECTORY( benchmark )
ADD_SUBDIRECTORY( sofft     )