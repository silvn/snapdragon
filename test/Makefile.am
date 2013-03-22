check_PROGRAMS = bvec_density serialize
AM_DEFAULT_SOURCE_EXT = .cpp
SOURCES = bvec_density.cpp serialize.cpp test.h
CC = g++
ACLOCAL_AMFLAGS = -I m4
AM_CPPFLAGS = -I../src -I$(top_srcdir)/../src -I$(top_srcdir)/../test
TESTS = $(check_PROGRAMS)
AM_LDFLAGS = -L$(top_srcdir)/src
if MACOS
    suff = -mt
endif
LDADD = -lsnapdragon \
	-lboost_serialization$(suff) \
	-lboost_thread$(suff) \
	-lboost_system$(suff)
