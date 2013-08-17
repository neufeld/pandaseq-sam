AC_INIT([pandaseq-sam], [1.0], [andre@masella.name])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_HEADERS(config.h)
AC_PROG_LIBTOOL
AC_HEADER_STDC

m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES([yes])], )

PLAIN_CPPFLAGS="$CPPFLAGS"
PLAIN_CFLAGS="$CFLAGS"
PLAIN_LIBS="$LIBS"

AC_CHECK_TYPES([bam1_t, khiter_t], [], [AC_MSG_ERROR([*** samtools is required, install the samtools sources in this directory. See the README for details])], [[
		#include "samtools/sam.h"
		#include "samtools/khash.h"
	]])

ACX_PTHREAD
AC_CHECK_HEADERS([math.h])
AC_CHECK_LIB([m],[pow])
AC_CHECK_HEADER([zlib.h], [], [AC_MSG_ERROR([*** zlib.h is required, install zlib header files])])
AC_CHECK_LIB([z], [gzopen], [], [AC_MSG_ERROR([*** gzopen is required, install zlib library files])])
LIB_CPPFLAGS="$CPPFLAGS"
LIB_CFLAGS="$CFLAGS"
LIB_LIBS="$LIBS"
AC_SUBST(LIB_CPPFLAGS)
AC_SUBST(LIB_CFLAGS)
AC_SUBST(LIB_LIBS)

PKG_CHECK_MODULES(PANDASEQ, [ pandaseq-2 >= 2.2 ])

CPPFLAGS="$PLAIN_CPPFLAGS"
CFLAGS="$PLAIN_CFLAGS"
LIBS="$PLAIN_LIBS"

LIB_NAME=pandaseq-sam-1
AC_SUBST(LIB_NAME)
# http://www.gnu.org/software/libtool/manual/html_node/Updating-version-info.html#Updating-version-info
LIB_VER=0:0:0
AC_SUBST(LIB_VER)
AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES([${LIB_NAME}.pc:${LIB_NAME}.pc.in], , [LIB_NAME=$LIB_NAME])
AC_CONFIG_FILES([build-macos-pkg:build-macos-pkg.in], [chmod +x build-macos-pkg])
AC_OUTPUT
