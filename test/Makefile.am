check_PROGRAMS = bvec_density serialize kmerize
AM_DEFAULT_SOURCE_EXT = .cpp
SOURCES = bvec_density.cpp serialize.cpp kmerize.cpp test.h
CC = g++
ACLOCAL_AMFLAGS = -I m4
AM_CXXFLAGS = -I../src -I$(top_srcdir)/../src -I$(top_srcdir)/../test
TESTS = $(check_PROGRAMS) histo.test
AM_LDFLAGS = -L$(top_srcdir)/src
if MACOS
    suff = -mt
endif
LDADD = -lsnapdragon \
	-lboost_serialization$(suff) \
	-lboost_thread$(suff) \
	-lboost_system$(suff)
TEST_LOG_DRIVER = env AM_TAP_AWK='$(AWK)' $(SHELL) \
	$(top_srcdir)/test/tap-driver.sh
EXTRA_DIST = $(TESTS)