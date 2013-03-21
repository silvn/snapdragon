bin_PROGRAMS = packmers kmerizer
packmers_SOURCES = packmers64.cpp
kmerizer_SOURCES = kmerizer.cpp
AM_CPPFLAGS = -I../countkmers
AM_LDFLAGS  = -L../countkmers/.libs
if MACOS
    suff=-mt
endif
LIBS = -lbvec -lz -lboost_thread$(suff) -lboost_system$(suff) \
	 -lboost_serialization$(suff)
