if MACOS
    suff = -mt
endif
LDADD = -lboost_serialization$(suff)
if BUILD_FASTBIT
   EXTRA_BUILD=fastbit
endif
AM_DEFAULT_SOURCE_EXT = .cpp
CC = g++

SUBDIRS = countkmers preprocessing kmerizer $(EXTRA_BUILD)
lib_LTLIBRARIES = libsnapdragon.la
libsnapdragon_la_SOURCES =
nodist_EXTRA_libsnapdragon_la_SOURCES = dummy.cxx
libsnapdragon_la_LIBADD = \
	countkmers/libbvec.la \
	kmerizer/libkmerizer.la
