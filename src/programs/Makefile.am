bin_PROGRAMS = histo kcount
histo_SOURCES = histo.cpp
kcount_SOURCES = kcount.cpp
AM_CXXFLAGS = -I../kmerizer -I../bvec
if MACOS
   suff = -mt
endif
AM_LDFLAGS = -L../kmerizer -L../bvec \
	-lkmerizer -lbvec -lz \
	-lboost_serialization$(suff) -lboost_system$(suff) -lboost_thread$(suff)
