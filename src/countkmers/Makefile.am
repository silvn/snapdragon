lib_LTLIBRARIES = libbvec.la
libbvec_la_SOURCES = bvec32.cpp bvec32.h kseq.h
AM_DEFAULT_SOURCE_EXT = .cpp
CC = g++
if MACOS
LDADD = -lboost_serialization-mt
else
LDADD = -lboost_serialization
endif
