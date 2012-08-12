AC_INIT([pandaseq-sam], [1.0], [andre@masella.name])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_HEADERS(config.h)
AC_PROG_LIBTOOL
AC_HEADER_STDC
AC_CHECK_HEADERS([math.h])
AC_CHECK_LIB([m],[pow])
ACX_PTHREAD

m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES([yes])], )

AC_CHECK_HEADER([zlib.h], [], [AC_MSG_ERROR([*** zlib.h is required, install zlib header files])])
AC_CHECK_LIB([z], [gzopen], [], [AC_MSG_ERROR([*** gzopen is required, install zlib library files])])
PKG_CHECK_MODULES(PANDASEQ, pandaseq-2)

LIB_NAME=pandaseq-sam-1
AC_SUBST(LIB_NAME)
AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES([${LIB_NAME}.pc:${LIB_NAME}.pc.in], , [LIB_NAME=$LIB_NAME])
AC_OUTPUT
