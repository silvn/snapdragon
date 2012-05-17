SUBDIRS = fastbit
bin_PROGRAMS = fbsql fbsql_tab join_once scatter gwas2go ranges stacker
SUFFIXES = .cpp
CC = g++
CPPFLAGS = -Ifastbit/src
LDFLAGS  = -Lfastbit/src/.libs -Wl,-rpath,$(libdir)
LIBS = -lfastbit -lm
if MACOS
stacker_LDADD = -lboost_thread-mt
else
stacker_LDADD = -lboost_thread
endif
