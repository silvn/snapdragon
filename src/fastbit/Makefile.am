bin_PROGRAMS = fbsql fbsql_tab join_once scatter gwas2go \
	ranges stacker loadBED playground
AM_DEFAULT_SOURCE_EXT = .cpp
CC = g++
AM_CPPFLAGS = -I$(srcdir)/vendor/fastbit/src
AM_LDFLAGS  = -L$(srcdir)/vendor/fastbit/src/.libs -Wl,-rpath,$(libdir)
LIBS = -lfastbit -lm
if MACOS
	suff=-mt
endif
stacker_LDADD = -lboost_thread$(suff) -lboost_system$(suff)
loadBED_LDADD = -lboost_thread$(suff) -lboost_system$(suff)
