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

# ==> Modules in this custom library [add the ones you have]:
set (SUPPL_PROGS_SRCS mbcontour.c mbswath.c mbgrdtiff.c mbclean_j.c mbgrid_j.c mbimport.c)
#=========================================================================
# Most likely no changes below here

# Find libraries
message (STATUS "Searching dependent libraries. This may take a few minutes...")

# lib sources
set (LIB_SRCS ${SUPPL_PROGS_SRCS})

# create build targets
string (REPLACE ".c" "" _gmt_modules "${SUPPL_PROGS_SRCS}")
install_module_symlink (${_gmt_modules})

# this prints out the final version of the compiler/linker flags
string (TOUPPER "${CMAKE_BUILD_TYPE}" CMAKE_BUILD_TYPE_UPPER)
string(REGEX REPLACE "[ ]+" " " _using_cflags ${CMAKE_C_FLAGS} " "
	${CMAKE_C_FLAGS_${CMAKE_BUILD_TYPE_UPPER}})
message (STATUS "Using CFLAGS = '${_using_cflags}'")

# vim: textwidth=78 noexpandtab tabstop=2 softtabstop=2 shiftwidth=2
