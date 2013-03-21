check_PROGRAMS = bvec_density serialize
AM_DEFAULT_SOURCE_EXT = .cpp
SOURCES = bvec_density.cpp serialize.cpp test.h
CC = g++
ACLOCAL_AMFLAGS = -I m4
AM_CPPFLAGS = -I../src -I$(top_srcdir)/../src -I$(top_srcdir)/../test
TESTS = $(check_PROGRAMS)
AM_LDFLAGS = -L../src/countkmers/.libs
LIBS = -lbvec
if MACOS
LDADD = -lboost_serialization-mt
else
LDADD = -lboost_serialization
endif
