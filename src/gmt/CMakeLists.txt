#
# $Id: CMakeLists.txt 11838 2013-06-28 06:49:32Z pwessel $
#
# Copyright (c) 1991-2015 by P. Wessel, W. H. F. Smith, R. Scharroo, J. Luis, and F. Wobbe
# See LICENSE.TXT file for copying and redistribution conditions.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation; version 3 or any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# Contact info: http://www.soest.hawaii.edu/PT/GSFML
#-------------------------------------------------------------------------------
#
# 1. To modify the cmake process: Edit your cmake/ConfigUser.cmake file
# 2. To add more modules: Edit CUSTOM_PROGS_SRCS below.
# 3. Edit this: LIB_STRING="GMT FZ: Tools for the GSFML project"

# This is src/CMakeLists.txt
#

# ==> Modules in this custom library [add the ones you have]:
set (CUSTOM_PROGS_SRCS mbcontour.c mbswath.c mbgrdtiff.c)
#=========================================================================
# Most likely no changes below here

include (ManageString)

# Include CREATE_DEBUG_SYM macro
include (CreateDebugSym)

# On Windows do not prefix the DLLs with the 'lib' prefix
if (WIN32)
	set(LIB_PREFIX "")
else (WIN32)
	set(LIB_PREFIX "lib")
endif (WIN32)

# Find libraries
message (STATUS "Searching dependent libraries. This may take a few minutes...")
find_package (GMT REQUIRED)
include_directories (${GMT_INCLUDE_DIR})
find_package (NETCDF REQUIRED)
include_directories (${NETCDF_INCLUDE_DIR})
include_directories (${ONCRPCINC})
include_directories (${GETOPTINC})
include_directories (${PROJINCDIR})

find_package (GDAL)
if (GDAL_FOUND)
	set (HAVE_GDAL TRUE CACHE INTERNAL "System has GDAL")
	include_directories (${GDAL_INCLUDE_DIR})
	list (APPEND GMT_OPTIONAL_LIBRARIES ${GDAL_LIBRARIES})
endif (GDAL_FOUND)

find_package (Threads)

# check for math and POSIX functions
include(ConfigureChecks)

if(NOT CUSTOM_SOURCE_DIR)
	set(CUSTOM_SOURCE_DIR "${CMAKE_SOURCE_DIR}")
endif()

include_directories(${CUSTOM_SOURCE_DIR}/include)
include_directories(${CMAKE_CURRENT_BINARY_DIR})
set(CMAKE_INCLUDE_CURRENT_DIR TRUE)

# Support code for the modules:
set (CUSTOM_LIB_SRCS ${CMAKE_PROJECT_NAME}_module.c)

message("\n MERDA  ${CMAKE_PROJECT_NAME}\n")

# lib targets
set (CUSTOM_LIBS mbgmtlib)

# mbgmtlib
set (CMAKE_SHARED_MODULE_SUFFIX .so)
add_library (mbgmtlib
	${CUSTOM_LIB_SRCS}
	${CUSTOM_PROGS_SRCS}
	)

target_link_libraries (mbgmtlib
	${GMT_LIBRARIES}
	${PSL_LIBRARY}
	${CMAKE_THREAD_LIBS_INIT}
	${PROJLIB}
	${ONCRPCLIB}
	${GETOPTLIB}
	${MBIOLIB}
	${MBAUXLIB}
	)

if (HAVE_M_LIBRARY)
	# link the math library
	target_link_libraries (mbgmtlib)
endif (HAVE_M_LIBRARY)

set_target_properties (mbgmtlib
	PROPERTIES
	OUTPUT_NAME ${CMAKE_PROJECT_NAME}
	RUNTIME_OUTPUT_NAME ${CMAKE_PROJECT_NAME}
	LIBRARY_OUTPUT_DIRECTORY plugins
	RUNTIME_OUTPUT_DIRECTORY plugins
	PREFIX ""
	LINK_FLAGS "${USER_CUSTOMLIB_LINK_FLAGS}"
	DEFINE_SYMBOL "LIBRARY_EXPORTS")

# If a renaming of the custom dll has been set in ConfigUser.cmake
if (WIN32 AND MBGMT_DLL_RENAME)
	set_target_properties (mbgmtlib PROPERTIES RUNTIME_OUTPUT_NAME ${MBGMT_DLL_RENAME})
endif ()

# generate Mac debugging symbols
#create_debug_sym (${CUSTOM_LIBS})

# temporary hack
if (NOT CUSTOM_INSTALL_MONOLITHIC)
	set (CUSTOM_MAIN_NAME ${CMAKE_PROJECT_NAME}${CUSTOM_INSTALL_NAME_SUFFIX})
else (NOT CUSTOM_INSTALL_MONOLITHIC)
	set (CUSTOM_MAIN_NAME ${CMAKE_PROJECT_NAME})
endif (NOT CUSTOM_INSTALL_MONOLITHIC)


install (CODE "
execute_process (
	COMMAND ${CMAKE_COMMAND} -E copy_if_different
	\"${CMAKE_CURRENT_SOURCE_DIR}/mbcontour\"
	\"\$ENV{DESTDIR}\${CMAKE_INSTALL_PREFIX}/${CUSTOM_BINDIR}/mbcontour${CUSTOM_INSTALL_NAME_SUFFIX}\")
execute_process (
	COMMAND ${CMAKE_COMMAND} -E copy_if_different
	\"${CMAKE_CURRENT_SOURCE_DIR}/mbswath\"
	\"\$ENV{DESTDIR}\${CMAKE_INSTALL_PREFIX}/${CUSTOM_BINDIR}/mbswath${CUSTOM_INSTALL_NAME_SUFFIX}\")
execute_process (
	COMMAND ${CMAKE_COMMAND} -E copy_if_different
	\"${CMAKE_CURRENT_SOURCE_DIR}/mbgrdtiff\"
	\"\$ENV{DESTDIR}\${CMAKE_INSTALL_PREFIX}/${CUSTOM_BINDIR}/mbgrdtiff${CUSTOM_INSTALL_NAME_SUFFIX}\")
")
# add the install targets
install (TARGETS ${CUSTOM_LIBS}
	ARCHIVE DESTINATION ${CUSTOM_LIBDIR}
	COMPONENT Runtime
	LIBRARY DESTINATION ${CUSTOM_LIBDIR}
	COMPONENT Runtime
	RUNTIME DESTINATION ${CUSTOM_BINDIR}
	COMPONENT Runtime)

if (WIN32 AND _debug_build)
	# install program database files
	set (_pdb_files)
	foreach (_target ${LIB_PREFIX}${CUSTOM_CUSTOM_RENAME})
		list (APPEND _pdb_files "${CMAKE_CURRENT_BINARY_DIR}/${_target}.pdb")
	endforeach ()
	install (FILES ${_pdb_files}
		DESTINATION ${CUSTOM_BINDIR}
		COMPONENT Runtime)
endif (WIN32 AND _debug_build)

# this prints out the final version of the compiler/linker flags
string (TOUPPER "${CMAKE_BUILD_TYPE}" CMAKE_BUILD_TYPE_UPPER)
string(REGEX REPLACE "[ ]+" " " _using_cflags ${CMAKE_C_FLAGS} " "
	${CMAKE_C_FLAGS_${CMAKE_BUILD_TYPE_UPPER}})
message (STATUS "Using CFLAGS = '${_using_cflags}'")

# vim: textwidth=78 noexpandtab tabstop=2 softtabstop=2 shiftwidth=2
