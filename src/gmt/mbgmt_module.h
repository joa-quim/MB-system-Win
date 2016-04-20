/*
 * $Id $
 * Copyright (c) 2015 by P. Wessel
 * See LICENSE.TXT file for copying and redistribution conditions.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation; version 3 or any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * Contact info: http://www.soest.hawaii.edu/PT/GSFML
 */

/* gsfml_module.h declares the prototypes for gsfml module functions
 * and the array that contains gsfml GMT module parameters such as name
 * and purpose strings.
 */

#pragma once
#ifndef _GMT_MBGMT_MODULE_H
#define _GMT_MBGMT_MODULE_H

#ifdef __cplusplus /* Basic C++ support */
extern "C" {
#endif

/* Declaration modifiers for DLL support (MSC et al) */
#include "declspec.h"

/* Prototypes of all modules in the GMT core library */
EXTERN_MSC int GMT_mbswath (void *API, int mode, void *args);
EXTERN_MSC int GMT_mbcontour (void *API, int mode, void *args);
EXTERN_MSC int GMT_mbgrdtiff (void *API, int mode, void *args);

/* Pretty print all modules in the GMT core library and their purposes */
EXTERN_MSC void gmt_mbgmt_module_show_all (void *API);

/* Undocumented API function for developers to get information about a module */
EXTERN_MSC const char * gmt_mbgmt_module_info (void *API, char *candidate);

#ifdef __cplusplus
}
#endif

#endif /* !_GMT_GSFML_MODULE_H */
