bin_PROGRAMS = fbsql fbsql_tab join_once scatter gwas2go \
	ranges stacker loadBED playground
AM_DEFAULT_SOURCE_EXT = .cpp
CC = g++
fastbit_dir = $(top_srcdir)/vendor/fastbit
AM_CPPFLAGS = -I$(fastbit_dir)/src
AM_LDFLAGS  = -L$(fastbit_dir)/src/.libs -Wl,-rpath,$(libdir)
LIBS = -lfastbit -lm
if MACOS
	suff=-mt
endif
stacker_LDADD = -lboost_thread$(suff) -lboost_system$(suff)
loadBED_LDADD = -lboost_thread$(suff) -lboost_system$(suff)
SUBDIRS = ../../vendor/fastbit #$(fastbit_dir)