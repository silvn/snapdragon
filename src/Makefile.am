SUBDIRS = fastbit
bin_PROGRAMS = fbsql fbsql_tab join_once scatter gwas2go ranges stacker
SUFFIXES = .cpp
CC = g++
CPPFLAGS = -Ifastbit/src
LDFLAGS  = -Lfastbit/src/.libs -Wl,-rpath,$(libdir)
LIBS = -lfastbit -lm -lboost_thread-mt
