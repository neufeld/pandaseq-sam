AC_INIT([pandaseq-sam], [1.2], [andre@masella.name])
AM_INIT_AUTOMAKE([-Wall -Werror foreign subdir-objects])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_HEADERS(config.h)
AC_ARG_ENABLE(threads, AC_HELP_STRING([--disable-threads], [disable thread support (default is autodetect)]))
AM_PROG_AR
AM_PROG_CC_C_O
AC_PROG_LIBTOOL
AC_HEADER_STDC

m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES([yes])], )

PLAIN_CPPFLAGS="$CPPFLAGS"
PLAIN_CFLAGS="$CFLAGS"
PLAIN_LIBS="$LIBS"

if test "$enable_threads" != no; then
	ACX_PTHREAD
else
	acx_pthread_ok=no
fi
AM_CONDITIONAL([PTHREAD], [test x$acx_pthread_ok = xyes])

# This is here because libhts does not correctly link again libm
AC_CHECK_LIB([m],[pow])

case "$host" in
*-mingw*)
	LIBS="$LIBS -lws2_32"
	NO_UNDEFINED='-no-undefined'
	CPPFLAGS="$CPPFLAGS -march=i486"
	;;
*)
	NO_UNDEFINED=''
esac
LIB_CPPFLAGS="$CPPFLAGS"
LIB_CFLAGS="$CFLAGS"
LIB_LIBS="$LIBS"
AC_SUBST(LIB_CPPFLAGS)
AC_SUBST(LIB_CFLAGS)
AC_SUBST(LIB_LIBS)
AC_SUBST(NO_UNDEFINED)

PKG_CHECK_MODULES(PANDASEQ, [ pandaseq-2 >= 2.5 ])
PKG_CHECK_MODULES(HTS, [ htslib ], [], [
	AC_CHECK_HEADER([htslib/sam.h], [], [AC_MSG_ERROR([*** htslib is required, install htslib header files])])
	AC_CHECK_LIB([hts], [hts_open], [], [AC_MSG_ERROR([*** htslib is required, install htslib library files])], "$PTHREAD_CFLAGS")
])

CPPFLAGS="$PLAIN_CPPFLAGS"
CFLAGS="$PLAIN_CFLAGS"
LIBS="$PLAIN_LIBS"

LIB_NAME=pandaseq-sam-1
AC_SUBST(LIB_NAME)
# http://www.gnu.org/software/libtool/manual/html_node/Updating-version-info.html#Updating-version-info
LIB_VER=1:0:0
AC_SUBST(LIB_VER)
AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES([${LIB_NAME}.pc:${LIB_NAME}.pc.in], , [LIB_NAME=$LIB_NAME])
AC_CONFIG_FILES([build-macos-pkg:build-macos-pkg.in], [chmod +x build-macos-pkg])
AC_OUTPUT
