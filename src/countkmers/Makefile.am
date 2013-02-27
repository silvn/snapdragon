lib_LTLIBRARIES = libbvec.la
libbvec_la_SOURCES = bvec32.cpp
AM_DEFAULT_SOURCE_EXT = .cpp
CC = g++
if MACOS
	LIBS = -lboost_serialization-mt
else
	LIBS = -lboost_serialization
endif
