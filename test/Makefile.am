TESTS = histo.test BitVectorSerialization BitVectorDensity Kmerizer
check_PROGRAMS = runtests
runtests_SOURCES = test-kmerizer.cpp test-bvec.cpp test.h
ACLOCAL_AMFLAGS = -I m4
GTEST_ROOT = $(top_srcdir)/vendor/googletest
AM_CXXFLAGS = -I$(top_srcdir)/../src -I$(GTEST_ROOT) -I$(GTEST_ROOT)/include
AM_LDFLAGS = -L$(top_srcdir)/src -L$(GTEST_ROOT)/lib
if MACOS
    suff = -mt
endif
LDADD = -lsnapdragon \
	-lboost_serialization$(suff) \
	-lboost_thread$(suff) \
	-lboost_system$(suff) \
	-lgtest
LOG_DRIVER = env AM_TAP_AWK='$(AWK)' $(SHELL) \
	$(top_srcdir)/test/tap-driver.sh
EXTRA_DIST = $(TESTS)