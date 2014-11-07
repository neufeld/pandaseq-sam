AC_INIT([pandaseq-sam], [1.3], [andre@masella.name])
AM_INIT_AUTOMAKE([-Wall -Werror foreign subdir-objects])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_HEADERS(config.h)
AC_ARG_ENABLE(threads, AC_HELP_STRING([--disable-threads], [disable thread support (default is autodetect)]))
AM_PROG_AR
AM_PROG_CC_C_O
AC_PROG_LIBTOOL
AC_HEADER_STDC

m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES([yes])], )

if test "$enable_threads" != no; then
	ACX_PTHREAD
else
	acx_pthread_ok=no
fi
AM_CONDITIONAL([PTHREAD], [test x$acx_pthread_ok = xyes])

PKG_CHECK_MODULES(PANDASEQ, [ pandaseq-2 >= 2.8 ])
PKG_CHECK_MODULES(HTS, [ htslib ], [], [
	ORIGINAL_CFLAGS="$CPPFLAGS"
	ORIGINAL_LIBS="$LIBS"
	# This is here because libhts does not correctly link against libm and pthread
	AC_CHECK_LIB([m],[pow])
	if test "$enable_threads" = no; then
		ACX_PTHREAD
	fi
	if test x$acx_pthread_ok = xno; then
		AC_MSG_ERROR([*** pthreads are needed by htslib])
	fi

	AC_CHECK_HEADER([htslib/sam.h], [], [AC_MSG_ERROR([*** htslib is required, install htslib header files])])
	AC_CHECK_LIB([hts], [hts_open], [], [AC_MSG_ERROR([*** htslib is required, install htslib library files])], [$PTHREAD_CFLAGS])
	HTS_CFLAGS="$CFLAGS $PTHREAD_CPPFLAGS"
	HTS_LIBS="$LIBS $PTHREAD_LIBS"
	AC_SUBST(HTS_CFLAGS)
	AC_SUBST(HTS_LIBS)
	CPPFLAGS="$ORIGINAL_CFLAGS"
	LIBS="$ORIGINAL_LIBS"
])

LIB_NAME=pandaseq-sam-1
AC_SUBST(LIB_NAME)
# http://www.gnu.org/software/libtool/manual/html_node/Updating-version-info.html#Updating-version-info
LIB_VER=1:0:0
AC_SUBST(LIB_VER)
AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES([${LIB_NAME}.pc:${LIB_NAME}.pc.in], , [LIB_NAME=$LIB_NAME])
AC_CONFIG_FILES([build-macos-pkg:build-macos-pkg.in], [chmod +x build-macos-pkg])
AC_OUTPUT
