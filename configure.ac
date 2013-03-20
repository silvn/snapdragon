AC_COPYRIGHT((c) 2013 Cold Spring Harbor Laboratory)

AC_INIT([snapdragon],[0.2],[olson@warelab.org])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])
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

AC_PROG_CC
AC_PROG_CXX
AC_PROG_LIBTOOL

AC_ARG_ENABLE(fastbit,
	[AS_HELP_STRING([--enable-fastbit],
		[Enable FastBit library and wrappers [default=no]])],
	[if test "$enableval" = "yes"; then fastbit="yes"; else fastbit="no"; fi],
	[fastbit="no"]
)
AC_MSG_RESULT(${fastbit})
AM_CONDITIONAL(BUILD_FASTBIT, test "$fastbit" = "yes")
AC_SUBST(BUILD_FASTBIT)

AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([
	Makefile
	src/Makefile
	src/fastbit/Makefile
])
if test "$fastbit" = "yes"; then
	AC_CONFIG_SUBDIRS([vendor/fastbit])
fi
AC_OUTPUT
AC_CACHE_SAVE

