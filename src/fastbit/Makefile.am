bin_PROGRAMS = fbsql fbsql_tab join_once scatter gwas2go \
	ranges stacker loadBED playground
AM_DEFAULT_SOURCE_EXT = .cpp
CC = g++
FASTBIT_DIR = ../../vendor/fastbit
AM_CPPFLAGS = -I$(FASTBIT_DIR)/src
AM_LDFLAGS  = -L$(FASTBIT_DIR)/src/.libs -Wl,-rpath,$(libdir)
LIBS = -lfastbit -lm
if MACOS
	suff=-mt
endif
stacker_LDADD = -lboost_thread$(suff) -lboost_system$(suff)
loadBED_LDADD = -lboost_thread$(suff) -lboost_system$(suff)
SUBDIRS = $(FASTBIT_DIR)