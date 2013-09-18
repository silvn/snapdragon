bin_PROGRAMS = kcount
kcount_SOURCES = kcount.cpp kmerizer.cpp
if MACOS
    suff=-mt
endif
AM_CPPFLAGS = -I../bvec

AM_LDFLAGS = -lz -lboost_serialization$(suff) -lboost_system$(suff) -lboost_thread$(suff)
