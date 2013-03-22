if MACOS
    suff=-mt
endif
AM_CPPFLAGS = -I../countkmers

bin_PROGRAMS = packmers
packmers_SOURCES = packmers64.cpp kmerizer.h
packmers_LDFLAGS  = -L../countkmers/.libs -lz -lboost_thread$(suff) \
	-lboost_system$(suff) -lboost_serialization$(suff)

noinst_LTLIBRARIES = libkmerizer.la
libkmerizer_la_SOURCES = kmerizer.cpp kmerizer.h
