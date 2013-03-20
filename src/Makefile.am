bin_PROGRAMS = fbsql fbsql_tab join_once scatter gwas2go ranges stacker loadBED playground
AM_DEFAULT_SOURCE_EXT = .cpp
CC = g++
AM_CPPFLAGS = -Ifastbit/src
AM_LDFLAGS  = -Lfastbit/src/.libs -Wl,-rpath,$(libdir)
LIBS = -lfastbit -lm
if MACOS
stacker_LDADD = -lboost_thread-mt -lboost_system-mt
loadBED_LDADD = -lboost_thread-mt -lboost_system-mt
else
stacker_LDADD = -lboost_thread -lboost_system
loadBED_LDADD = -lboost_thread -lboost_system
endif
