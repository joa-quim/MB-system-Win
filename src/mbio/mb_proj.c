/*--------------------------------------------------------------------
 *    The MB-system:	mb_proj.c	7/16/2002
 *    $Id: mb_proj.c 2308 2017-06-04 19:55:48Z caress $
 *
 *    Copyright (c) 2002-2017 by
 *    David W. Caress (caress@mbari.org)
 *      Monterey Bay Aquarium Research Institute
 *      Moss Landing, CA 95039
 *    and Dale N. Chayes (dale@ldeo.columbia.edu)
 *      Lamont-Doherty Earth Observatory
 *      Palisades, NY 10964
 *
 *    See README file for copying and redistribution conditions.
 *--------------------------------------------------------------------*/
/*
 * mb_proj.c includes the "mb_" functions used to initialize
 * projections, and then to do forward (mb_proj_forward())
 * and inverse (mb_proj_inverse()) projections
 * between geographic coordinates (longitude and latitude) and
 * projected coordinates (e.g. eastings and northings in meters).
 * One can also tranlate between coordinate systems using mb_proj_transform().
 * This code uses libproj. The code in libproj derives without modification
 * from the PROJ.4 distribution. PROJ was originally developed by
 * Gerard Evandim, and is now maintained and distributed by
 * Frank Warmerdam, <warmerdam@pobox.com>
 *
 * David W. Caress
 * July 16, 2002
 * RVIB Nathaniel B. Palmer
 * Somewhere west of Conception, Chile
 *
 * Author:	D. W. Caress
 * Date:	July 16, 2002
 *
 *
 */

/* standard include files */
#include <stdio.h>
#include <math.h>
#include <unistd.h>
#include <sys/stat.h>

/* mbio include files */
#include "mb_status.h"
#include "mb_define.h"
#include "proj_api.h"

#ifndef _WIN32
#include "projections.h"
#else
char *GMT_runtime_bindir_win32(char *result);
#endif

static char rcs_id[] = "$Id: mb_proj.c 2308 2017-06-04 19:55:48Z caress $";

#ifdef _WIN32
char *GMT_runtime_bindir_win32 (char *result);
#endif

/*--------------------------------------------------------------------*/
int mb_proj_init(int verbose, char *projection, void **pjptr, int *error) {
	char *function_name = "mb_proj_init";
	int status = MB_SUCCESS;
	char pj_init_args[MB_PATH_MAXLINE];
	projPJ pj;
	struct stat file_status;
	int	fstat;
#ifndef _WIN32
#	include "projections.h"
#else
	/* But on Windows get it from the bin dir */
#	include <unistd.h>
	char *pch, projectionfile[MB_PATH_MAXLINE + 1];
	/* Find the path to the bin directory and from it, the location of the Projections.dat file */
	GMT_runtime_bindir_win32 (projectionfile);
	pch = strrchr(projectionfile, '\\');		/* Seek for the last '\' or '/'. One of them must exist. */
	if (pch == NULL)
		pch = strrchr(projectionfile, '/');
	pch[0] = '\0';
	strcat(projectionfile, "\\share\\mbsystem\\Projections.dat");
#endif

	/* print input debug statements */
	if (verbose >= 2) {
		fprintf(stderr, "\ndbg2  MBIO function <%s> called\n", function_name);
		fprintf(stderr, "dbg2  Revision id: %s\n", rcs_id);
		fprintf(stderr, "dbg2  Input arguments:\n");
		fprintf(stderr, "dbg2       verbose:    %d\n", verbose);
		fprintf(stderr, "dbg2       projection: %s\n", projection);
	}

/* Normally the header file projections.h sets the location of the
	projections.dat file in a string projectionfile, but on Windows instead
	use GMT constructs to find the path to the bin directory and from it,
	the location of the Projections.dat file.
	This construct has been defined by Joaquim Luis. */
#ifdef _WIN32
	GMT_runtime_bindir_win32(projectionfile);
	pch = strrchr(projectionfile, '\\'); /* Seek for the last '\' or '/'. One of them must exist. */
	if (pch == NULL)
		pch = strrchr(projectionfile, '/');
	pch[0] = '\0';
	strcat(projectionfile, "\\share\\mbsystem\\Projections.dat");
#endif

	/* check the existence of the projection database */
	if ((fstat = stat(projectionfile, &file_status)) == 0 && (file_status.st_mode & S_IFMT) != S_IFDIR) {
		/* initialize the projection */
		sprintf(pj_init_args, "+init=%s:%s", projectionfile, projection);
		pj = pj_init_plus(pj_init_args);
		*pjptr = (void *)pj;

		/* check success */
		if (*pjptr != NULL) {
			*error = MB_ERROR_NO_ERROR;
			status = MB_SUCCESS;
		}
		else {
			*error = MB_ERROR_BAD_PROJECTION;
			status = MB_FAILURE;
		}
	}
	else {
		/* cannot initialize the projection */
		fprintf(stderr, "\nUnable to open projection database at expected location:\n\t%s\n", projectionfile);
		fprintf(stderr, "Set projection database location using the $MBSYSTEM_HOME and $PROJECTIONS \ntags in the "
		                "install_makefiles script.\n\n");
		*error = MB_ERROR_MISSING_PROJECTIONS;
		status = MB_FAILURE;
	}

	/* print output debug statements */
	if (verbose >= 2) {
		fprintf(stderr, "\ndbg2  MBIO function <%s> completed\n", function_name);
		fprintf(stderr, "dbg2  Revision id: %s\n", rcs_id);
		fprintf(stderr, "dbg2  Return values:\n");
		fprintf(stderr, "dbg2       pjptr:           %p\n", (void *)*pjptr);
		fprintf(stderr, "dbg2       error:           %d\n", *error);
		fprintf(stderr, "dbg2  Return status:\n");
		fprintf(stderr, "dbg2       status:          %d\n", status);
	}

	/* return status */
	return (status);
}
/*--------------------------------------------------------------------*/
int mb_proj_free(int verbose, void **pjptr, int *error) {
	char *function_name = "mb_proj_free";
	int status = MB_SUCCESS;
	projPJ pj;

	/* print input debug statements */
	if (verbose >= 2) {
		fprintf(stderr, "\ndbg2  MBIO function <%s> called\n", function_name);
		fprintf(stderr, "dbg2  Revision id: %s\n", rcs_id);
		fprintf(stderr, "dbg2  Input arguments:\n");
		fprintf(stderr, "dbg2       verbose:    %d\n", verbose);
		fprintf(stderr, "dbg2       pjptr:      %p\n", (void *)*pjptr);
	}

	/* free the projection */
	if (pjptr != NULL) {
		pj = (projPJ)*pjptr;
		pj_free(pj);
		*pjptr = NULL;
	}

	/* assume success */
	*error = MB_ERROR_NO_ERROR;
	status = MB_SUCCESS;

	/* print output debug statements */
	if (verbose >= 2) {
		fprintf(stderr, "\ndbg2  MBIO function <%s> completed\n", function_name);
		fprintf(stderr, "dbg2  Revision id: %s\n", rcs_id);
		fprintf(stderr, "dbg2  Return values:\n");
		fprintf(stderr, "dbg2       pjptr:           %p\n", (void *)*pjptr);
		fprintf(stderr, "dbg2       error:           %d\n", *error);
		fprintf(stderr, "dbg2  Return status:\n");
		fprintf(stderr, "dbg2       status:          %d\n", status);
	}

	/* return status */
	return (status);
}
/*--------------------------------------------------------------------*/
int mb_proj_forward(int verbose, void *pjptr, double lon, double lat, double *easting, double *northing, int *error) {
	char *function_name = "mb_proj_forward";
	int status = MB_SUCCESS;
	projPJ pj;
	projUV pjxy;
	projUV pjll;

	/* print input debug statements */
	if (verbose >= 2) {
		fprintf(stderr, "\ndbg2  MBIO function <%s> called\n", function_name);
		fprintf(stderr, "dbg2  Revision id: %s\n", rcs_id);
		fprintf(stderr, "dbg2  Input arguments:\n");
		fprintf(stderr, "dbg2       verbose:    %d\n", verbose);
		fprintf(stderr, "dbg2       pjptr:      %p\n", (void *)pjptr);
		fprintf(stderr, "dbg2       lon:        %f\n", lon);
		fprintf(stderr, "dbg2       lat:        %f\n", lat);
	}

	/* do forward projection */
	if (pjptr != NULL) {
		pj = (projPJ)pjptr;
		pjll.u = DTR * lon;
		pjll.v = DTR * lat;
		pjxy = pj_fwd(pjll, pj);
		*easting = pjxy.u;
		*northing = pjxy.v;
	}

	/* assume success */
	*error = MB_ERROR_NO_ERROR;
	status = MB_SUCCESS;

	/* print output debug statements */
	if (verbose >= 2) {
		fprintf(stderr, "\ndbg2  MBIO function <%s> completed\n", function_name);
		fprintf(stderr, "dbg2  Revision id: %s\n", rcs_id);
		fprintf(stderr, "dbg2  Return values:\n");
		fprintf(stderr, "dbg2       easting:         %f\n", *easting);
		fprintf(stderr, "dbg2       northing:        %f\n", *northing);
		fprintf(stderr, "dbg2       error:           %d\n", *error);
		fprintf(stderr, "dbg2  Return status:\n");
		fprintf(stderr, "dbg2       status:          %d\n", status);
	}

	/* return status */
	return (status);
}
/*--------------------------------------------------------------------*/
int mb_proj_inverse(int verbose, void *pjptr, double easting, double northing, double *lon, double *lat, int *error) {
	char *function_name = "mb_proj_inverse";
	int status = MB_SUCCESS;
	projPJ pj;
	projUV pjxy;
	projUV pjll;

	/* print input debug statements */
	if (verbose >= 2) {
		fprintf(stderr, "\ndbg2  MBIO function <%s> called\n", function_name);
		fprintf(stderr, "dbg2  Revision id: %s\n", rcs_id);
		fprintf(stderr, "dbg2  Input arguments:\n");
		fprintf(stderr, "dbg2       verbose:    %d\n", verbose);
		fprintf(stderr, "dbg2       pjptr:      %p\n", (void *)pjptr);
		fprintf(stderr, "dbg2       easting:    %f\n", easting);
		fprintf(stderr, "dbg2       northing:   %f\n", northing);
	}

	/* do forward projection */
	if (pjptr != NULL) {
		pj = (projPJ)pjptr;
		pjxy.u = easting;
		pjxy.v = northing;
		pjll = pj_inv(pjxy, pj);
		*lon = RTD * pjll.u;
		*lat = RTD * pjll.v;
	}

	/* assume success */
	*error = MB_ERROR_NO_ERROR;
	status = MB_SUCCESS;

	/* print output debug statements */
	if (verbose >= 2) {
		fprintf(stderr, "\ndbg2  MBIO function <%s> completed\n", function_name);
		fprintf(stderr, "dbg2  Revision id: %s\n", rcs_id);
		fprintf(stderr, "dbg2  Return values:\n");
		fprintf(stderr, "dbg2       lon:             %f\n", *lon);
		fprintf(stderr, "dbg2       lat:             %f\n", *lat);
		fprintf(stderr, "dbg2       error:           %d\n", *error);
		fprintf(stderr, "dbg2  Return status:\n");
		fprintf(stderr, "dbg2       status:          %d\n", status);
	}

	/* return status */
	return (status);
}
/*--------------------------------------------------------------------*/
int mb_proj_transform(int verbose, void *pjsrcptr, void *pjdstptr, int npoint, double *x, double *y, double *z, int *error) {
	char *function_name = "mb_proj_transform";
	int status = MB_SUCCESS;
	int i;

	/* print input debug statements */
	if (verbose >= 2) {
		fprintf(stderr, "\ndbg2  MBIO function <%s> called\n", function_name);
		fprintf(stderr, "dbg2  Revision id: %s\n", rcs_id);
		fprintf(stderr, "dbg2  Input arguments:\n");
		fprintf(stderr, "dbg2       verbose:    %d\n", verbose);
		fprintf(stderr, "dbg2       pjptr:      %p\n", (void *)pjsrcptr);
		fprintf(stderr, "dbg2       pjptr:      %p\n", (void *)pjdstptr);
		fprintf(stderr, "dbg2       npoint:     %d\n", npoint);
		for (i = 0; i < npoint; i++)
			fprintf(stderr, "dbg2       point[%d]:  x:%f y:%f z:%f\n", i, x[i], y[i], z[i]);
	}

	/* do transform */
	if (pjsrcptr != NULL && pjdstptr != NULL) {
		pj_transform((projPJ *)pjsrcptr, (projPJ *)pjdstptr, npoint, 1, x, y, z);
	}

	/* assume success */
	*error = MB_ERROR_NO_ERROR;
	status = MB_SUCCESS;

	/* print output debug statements */
	if (verbose >= 2) {
		fprintf(stderr, "\ndbg2  MBIO function <%s> completed\n", function_name);
		fprintf(stderr, "dbg2  Revision id: %s\n", rcs_id);
		fprintf(stderr, "dbg2  Return values:\n");
		fprintf(stderr, "dbg2       npoint:     %d\n", npoint);
		for (i = 0; i < npoint; i++)
			fprintf(stderr, "dbg2       point[%d]:  x:%f y:%f z:%f\n", i, x[i], y[i], z[i]);
		fprintf(stderr, "dbg2       error:           %d\n", *error);
		fprintf(stderr, "dbg2  Return status:\n");
		fprintf(stderr, "dbg2       status:          %d\n", status);
	}

	/* return status */
	return (status);
}
/*--------------------------------------------------------------------*/
