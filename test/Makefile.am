GTEST_DIR = $(top_srcdir)/vendor/googletest

check_PROGRAMS = test-kmerizer test-bvec test-freqmap
TESTS = histo.test $(check_PROGRAMS)
AM_DEFAULT_SOURCE_EXT = .cpp
CC = g++
SOURCES = test-kmerizer.cpp test-bvec.cpp test-freqmap.cpp test.h
ACLOCAL_AMFLAGS = -I m4

all-local: $(GTEST_DIR)/make/libgtest.a

$(GTEST_DIR)/make/libgtest.a:
	@ echo "Compiling GoogleTest"
	@ cd $(GTEST_DIR)/make; \
		$(MAKE) $(AM_MAKEFLAGS) all; \
		ar -rv libgtest.a gtest-all.o

clean-local:
	@ echo "Making clean in GoogleTest";
	@ rm libgtest.a
	@ cd $(GTEST_DIR)/make; \
		$(MAKE) $(AM_MAKEFLAGS) clean


AM_CXXFLAGS = -I$(top_srcdir)/src -I$(GTEST_DIR)/include
AM_LDFLAGS = -L$(top_srcdir)/src -L$(GTEST_DIR)/make -isystem $(GTEST_DIR)/include -pthread
if MACOS
    suff = -mt
endif
LDADD = -lsnapdragon \
	-lboost_serialization$(suff) \
	-lboost_thread$(suff) \
	-lboost_system$(suff) \
	-lgtest

EXTRA_DIST = $(TESTS)
