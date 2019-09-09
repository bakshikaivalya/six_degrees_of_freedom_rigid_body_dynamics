# CMAKE file for the mesh directory.

CMAKE_MINIMUM_REQUIRED (VERSION 2.6)
#SET(CMAKE_VERBOSE_MAKEFILE On) 

# Set the list of sources for this library.    
SET( sources
    define_variables.cpp
    calculate_structural_prop.cpp
)

# include directories
INCLUDE_DIRECTORIES( . )

# compiler flags    
SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pedantic -Weffc++ " )
SET( CMAKE_CXX_FLAGS_DEBUG " -g -pg" )
SET( CMAKE_CXX_FLAGS_RELEASE " -s -O3" )

# The target library to be built.
ADD_EXECUTABLE( structTool ${sources} )
