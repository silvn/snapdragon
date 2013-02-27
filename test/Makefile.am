check_PROGRAMS = bvec_density serialize
AM_DEFAULT_SOURCE_EXT = .cpp
CC = g++
ACLOCAL_AMFLAGS = -I m4
AM_CPPFLAGS = -I../src
TESTS = $(check_PROGRAMS)
AM_LDFLAGS = -L../src/countkmers/.libs
LIBS = -lbvec
if MACOS
LDADD = -lboost_serialization-mt
else
LDADD = -lboost_serialization
endif
