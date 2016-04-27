/*
 * $Id: mbgmt_module.c 352 2015-02-01 04:34:58Z pwessel $
 *
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
 *--------------------------------------------------------------------
 *
 * mbgmt_module.c populates the local array g_mbgmt_module
 * with parameters such as name, group and purpose strings.
 * This file also contains the following convenience function to
 * display all module purposes:
 *
 *   void gmt_mbgmt_module_show_all (struct GMTAPI_CTRL *API);
 *
 * This function will be called by gmt --help
 */

#include "gmt.h"
#include "mbgmt_module.h"
#include <string.h>

/* Sorted array with information for all mbgmt modules */

/* name, library, and purpose for each module */
struct Gmt_moduleinfo {
	const char *name;             /* Program name */
	const char *component;        /* Component (core, supplement, custom) */
	const char *purpose;          /* Program purpose */
	const char *keys;             /* Program option info for external APIs */
};

struct Gmt_moduleinfo g_mbgmt_module[] = {
	{"mbcontour", "mbgmt", "Plot contours", ""},
	{"mbswath",   "mbgmt", "Plot swath bathymetry, amplitude, or backscatter", ""},
	{"mbgrdtiff", "mbgmt", "Create geotiff images", ""},
	{"mbclean_j", "mbgmt", "Identifies and flags artifacts in swath sonar bathymetry data", ""},
	{"mbgrid_j",  "mbgmt", "grid bathymetry, amplitude, or sidescan data of a set of swath sonar data files", ""},
	{"mbimport",  "mbgmt", "Get swath bathymetry, amplitude, or backscatter Image into Matlab", "<D{,CC(,MI}"},
	{NULL, NULL, NULL} /* last element == NULL detects end of array */
};

/* Pretty print all GMT mbgmt module names and their purposes */
void gmt_mbgmt_module_show_all (void *V_API) {
	unsigned int module_id = 0;
	char message[256];
	GMT_Message (V_API, GMT_TIME_NONE, "\n=== " "GMT mbgmt : Tools for the MBGMT project" " ===\n");
	while (g_mbgmt_module[module_id].name != NULL) {
		if (module_id == 0 || strcmp (g_mbgmt_module[module_id-1].component, g_mbgmt_module[module_id].component)) {
			/* Start of new supplemental group */
			sprintf (message, "\nModule name:     Purpose of %s module:\n", g_mbgmt_module[module_id].component);
			GMT_Message (V_API, GMT_TIME_NONE, message);
			GMT_Message (V_API, GMT_TIME_NONE, "----------------------------------------------------------------\n");
		}
		sprintf(message, "%-16s %s\n", g_mbgmt_module[module_id].name, g_mbgmt_module[module_id].purpose);
		GMT_Message (V_API, GMT_TIME_NONE, message);
		module_id++;
	}
}

/* Lookup module id by name, return option keys pointer (for external API developers) */
const char *gmt_mbgmt_module_info (void *API, char *candidate) {
	int module_id = 0;

	/* Match actual_name against g_module[module_id].name */
	while (g_mbgmt_module[module_id].name != NULL && strcmp (candidate, g_mbgmt_module[module_id].name))
		module_id++;

	/* Return Module keys or NULL */
	return (g_mbgmt_module[module_id].keys);
}
