#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([bvec], [0.0.1], [aolson@me.com])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])
AC_REVISION([m4_esyscmd_s([git describe --always])])

AC_CONFIG_SRCDIR([config.h.in])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_MACRO_DIR([m4])
AM_CONDITIONAL([MACOS], [test "$(uname -s)" = "Darwin"])
# Options to configure.
# Boost
AC_ARG_WITH(boost, AS_HELP_STRING([--with-boost=PATH],
    [directory for the boost header files]))
if test "$with_boost" -a -d "$with_boost"; then
    boost_cppflags="-I$with_boost"
fi
AC_SUBST(CPPFLAGS, "$boost_cppflags $CPPFLAGS")

AC_LANG_PUSH([C++])
AC_CHECK_HEADERS([boost/foreach.hpp], [],
    [AC_MSG_ERROR(You need the Boost libraries.)])
AC_LANG_POP([C++])

LIB_VERSION_INFO=0:9:0

AM_PROG_AR

# Checks for programs.
AC_PROG_CXX
AC_PROG_AWK
AC_PROG_CC
AC_PROG_CPP
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MAKE_SET
AC_PROG_LIBTOOL

# Checks for libraries.
AC_CHECK_LIB([z], [deflateInit])
AC_CHECK_LIB([boost_serialization], [serialize])

# Checks for header files.
AC_CHECK_HEADERS([stdlib.h string.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_SIZE_T
AC_TYPE_UINT32_T
AC_TYPE_UINT64_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([mkdir strerror])

AC_CONFIG_FILES([Makefile
                 src/Makefile
                 src/countkmers/Makefile
                 src/preprocessing/Makefile
                 test/Makefile])
LT_INIT
AC_OUTPUT
