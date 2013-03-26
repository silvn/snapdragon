bin_PROGRAMS = histo
histo_SOURCES = histo.cpp
AM_CXXFLAGS = -I../kmerizer -I../countkmers
if MACOS
   suff = -mt
endif
AM_LDFLAGS = -L../kmerizer -L../countkmers \
	-lkmerizer -lbvec -lz \
	-lboost_serialization$(suff) -lboost_system$(suff) -lboost_thread$(suff)
