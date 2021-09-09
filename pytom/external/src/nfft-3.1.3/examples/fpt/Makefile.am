# $Id: Makefile.am 3110 2009-03-13 16:32:18Z keiner $

AM_CPPFLAGS = -I$(top_srcdir)/include

exampledir = $(pkgdatadir)/examples/fpt

# programs
example_PROGRAMS = simple_test

# sources
simple_test_SOURCES = simple_test.c
simple_test_LDADD = $(top_builddir)/libnfft3.la
