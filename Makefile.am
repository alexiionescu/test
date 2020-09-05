ACLOCAL_AMFLAGS = -I m4

lib_LTLIBRARIES =
noinst_LTLIBRARIES =
bin_PROGRAMS =
noinst_PROGRAMS =
dist_bin_SCRIPTS =
CLEANFILES =
EXTRA_DIST =
pkginclude_HEADERS =

lib_LTLIBRARIES += lib/libgenetics.la
pkginclude_HEADERS += lib/genetics/genetics.h
lib_libgenetics_la_SOURCES = lib/genetics/genetics.h lib/genetics/genetics.c \
                       lib/genetics/transl_table.h lib/genetics/transl_table.c 

bin_PROGRAMS += bin/testam
bin_testam_SOURCES = src/main.c \
        src/test.h src/tests.c
bin_testam_LDADD = lib/libgenetics.la

dist_doc_DATA = README
