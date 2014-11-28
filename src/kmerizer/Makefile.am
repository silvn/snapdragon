bin_PROGRAMS = kcount
kcount_SOURCES = kcount.cpp
noinst_LTLIBRARIES = libkmerizer.la
libkmerizer_la_SOURCES = kmerizer.tpp kmerizer.h packmers64.cpp
AM_DEFAULT_SOURCE_EXT = .tpp
CC = g++
if MACOS
    suff=-mt
endif
AM_CPPFLAGS = -I../bitmap

AM_LDFLAGS = -lz -lboost_system$(suff) -lboost_thread$(suff)
