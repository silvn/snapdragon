AC_COPYRIGHT((c) 2012 Cold Spring Harbor Laboratory)

AC_INIT([snapdragon],[1.0],[olson@warelab.org])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])
AC_LANG_PUSH([C++])
AM_CONDITIONAL([MACOS], [test "$(uname -s)" = "Darwin"])
AC_LANG_POP([C++])
# Options to configure.
# Boost
AC_ARG_WITH(boost, AS_HELP_STRING([--with-boost=PATH],
    [directory for the boost header files]))
if test "$with_boost" -a -d "$with_boost"; then
    boost_cppflags="-I$with_boost"
fi
AC_SUBST(CPPFLAGS, "$boost_cppflags $CPPFLAGS")
AC_CHECK_HEADERS([boost/foreach.hpp], [],
    [AC_MSG_ERROR(You need the Boost libraries.)])
AC_PROG_CC
AC_PROG_CXX
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([
	Makefile
	src/Makefile
])
AC_CONFIG_SUBDIRS([src/fastbit])
AC_OUTPUT
