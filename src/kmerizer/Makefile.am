bin_PROGRAMS = kcount
kcount_SOURCES = kcount.cpp
if MACOS
    suff=-mt
endif
AM_CPPFLAGS = -I../bitmap

# AM_LDFLAGS = -lz -lboost_serialization$(suff) -lboost_system$(suff) -lboost_thread$(suff)
AM_LDFLAGS = -lz -lboost_system$(suff) -lboost_thread$(suff)
