/* src/mbio/mb_config.h.in.  Generated from configure.ac by autoheader.  */

/* Machine is littleendian, (Byteswapping on) */
#define BYTESWAPPED 1

/* Machine is bigendian, (Byteswapping off) */
/* #undef ENDIAN_BIG */

/* Turned on OpenGL define in config */
#undef GOT_GL

/* Turned on Motif define in config */
//#undef GOT_MOTIF

/* Define to 1 if you have the <dlfcn.h> header file. */
#undef HAVE_DLFCN_H

/* Define to 1 if you have the <GL/glu.h> header file. */
#undef HAVE_GL_GLU_H

/* Define to 1 if you have the <GL/glx.h> header file. */
#undef HAVE_GL_GLX_H

/* Define to 1 if you have the <GL/gl.h> header file. */
#undef HAVE_GL_GL_H

/* Define to 1 if you have the <inttypes.h> header file. */
//#undef HAVE_INTTYPES_H

/* Define to 1 if you have the `m' library (-lm). */
#undef HAVE_LIBM

/* Define to 1 if you have the `ws2_32' library (-lws2_32). */
#undef HAVE_LIBWS2_32

/* Have malloc.h */
//#undef HAVE_MALLOC_H

/* Define to 1 if you have the <memory.h> header file. */
//#undef HAVE_MEMORY_H

/* Define to 1 if you have the <rpc/rpc.h> header file. */
//#undef HAVE_RPC_RPC_H
#define HAVE_RPC_RPC_H 1

/* Define to 1 if you have the <rpc/types.h> header file. */
//#undef HAVE_RPC_TYPES_H
#define HAVE_RPC_TYPES_H 1

/* Define to 1 if you have the <stdint.h> header file. */
//#undef HAVE_STDINT_H

/* Define to 1 if you have the <stdlib.h> header file. */
//#undef HAVE_STDLIB_H

/* Define to 1 if you have the <strings.h> header file. */
//#undef HAVE_STRINGS_H

/* Define to 1 if you have the <string.h> header file. */
#undef HAVE_STRING_H

/* Define to 1 if you have the <sys/stat.h> header file. */
#undef HAVE_SYS_STAT_H

/* Define to 1 if you have the <sys/types.h> header file. */
#undef HAVE_SYS_TYPES_H

/* Define to 1 if you have the <unistd.h> header file. */
#undef HAVE_UNISTD_H

/* Define to the sub-directory in which libtool stores uninstalled libraries.
   */
#undef LT_OBJDIR

/* Define to 1 if your C compiler doesn't accept -c and -o together. */
#undef NO_MINUS_C_MINUS_O

/* Build libmbtrn and mbtrnpreprocess */
//#define MBTRN_ENABLED 1

/* Name of package */
#undef PACKAGE

/* Define to the address where bug reports for this package should be sent. */
#undef PACKAGE_BUGREPORT

/* Define to the full name of this package. */
#undef PACKAGE_NAME

/* Define to the full name and version of this package. */
#define PACKAGE_STRING "mbsystem 5.5.2350"

/* Define to the one symbol short name of this package. */
#undef PACKAGE_TARNAME

/* Define to the home page for this package. */
#undef PACKAGE_URL

/* Define to the version of this package. */
#define PACKAGE_VERSION "5.5.2350"

/* Define to 1 if you have the ANSI C header files. */
#undef STDC_HEADERS

/* Version number of package */
#define VERSION "5.5.2350"

/* Set VERSION_DATE define in mb_config.h */
#define VERSION_DATE "6 September 2018"

/* Turned on WIN32 define in config */
//#undef WIN32

/* Build with GSF */
//#undef WITH_GSF
#define WITH_GSF

/* Define WORDS_BIGENDIAN to 1 if your processor stores words with the most
   significant byte first (like Motorola and SPARC, unlike Intel). */
#if defined AC_APPLE_UNIVERSAL_BUILD
# if defined __BIG_ENDIAN__
#  define WORDS_BIGENDIAN 1
# endif
#else
# ifndef WORDS_BIGENDIAN
#  undef WORDS_BIGENDIAN
# endif
#endif

/* Define to the type of a signed integer type of width exactly 8 bits if such
   a type exists and the standard includes do not define it. */
//#undef int8_t
