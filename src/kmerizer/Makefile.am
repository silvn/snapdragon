bin_PROGRAMS = packmers
packmers_SOURCES = packmers64.cpp
if MACOS
    suff=-mt
endif
LIBS = -lz -lboost_thread$(suff) -lboost_system$(suff)
