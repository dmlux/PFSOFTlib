# - Find FFTW
# Find the native FFTW includes and library
#
#  FFTW_INCLUDE_DIR - where to find fftw3.h
#  FFTW_LIB   		- List of libraries when using FFTW.
#  FFTW_FOUND       - True if FFTW found.

IF(FFTW_INCLUDES)
    # Already in cache, be silent
    SET(FFTW_FIND_QUIETLY TRUE)
ENDIF(FFTW_INCLUDES)

SET(FFTW_INCLUDE_SEARCH_PATHS
  	/usr/include
  	/usr/local/include
  	$ENV{FFTW_HOME}
  	$ENV{FFTW_HOME}/include
)

SET(FFTW_LIB_SEARCH_PATHS
    /lib64/
    /lib/
    /usr/lib64
    /usr/lib
	/usr/local/lib64
	/usr/local/lib
	/usr/local/Cellar/fftw/3.3.4_1/lib/
	$ENV{FFTW}
	$ENV{FFTW}/lib
	$ENV{FFTW_HOME}
	$ENV{FFTW_HOME}/lib
)

FIND_PATH(FFTW_INCLUDE_DIR fftw3.h PATHS ${FFTW_INCLUDE_SEARCH_PATHS} SHARED IMPORTED)
FIND_LIBRARY(FFTW_LIB NAMES fftw3 PATHS ${FFTW_LIB_SEARCH_PATHS})

# handle the QUIETLY and REQUIRED arguments and set FFTW_FOUND to TRUE if
# all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(FFTW DEFAULT_MSG FFTW_LIB FFTW_INCLUDE_DIR)

# Check include files
IF(NOT FFTW_INCLUDE_DIR)
	SET(FFTW_FOUND OFF)
    MESSAGE(STATUS "~> Could not find FFTW include. Turning FFTW_FOUND off")
ENDIF()

# Check libraries
IF(NOT FFTW_LIB)
    SET(FFTW_FOUND OFF)
    MESSAGE(STATUS "~> Could not find FFTW lib. Turning FFTW_FOUND off")
ENDIF()

IF (FFTW_FOUND)
	IF (NOT FFTW_FIND_QUIETLY)
    	MESSAGE(STATUS "~> Found FFTW libraries: ${FFTW_LIB}")
		MESSAGE(STATUS "~> Found FFTW include: ${FFTW_INCLUDE_DIR}")
	ENDIF (NOT FFTW_FIND_QUIETLY)
ELSE (FFTW_FOUND)
	IF (FFTW_FIND_REQUIRED)
    	MESSAGE(FATAL_ERROR "~> Could not find FFTW")
	ENDIF (FFTW_FIND_REQUIRED)
ENDIF (FFTW_FOUND)

MARK_AS_ADVANCED(
	FFTW_LIBRARIES 
	FFTW_INCLUDES
)