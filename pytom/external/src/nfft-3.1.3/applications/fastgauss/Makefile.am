# $Id: Makefile.am 3086 2009-03-11 07:36:56Z keiner $

AM_CPPFLAGS = -I$(top_srcdir)/include

applicationdir = $(pkgdatadir)/applications/fastgauss

application_PROGRAMS = fastgauss

fastgauss_SOURCES = fastgauss.c
fastgauss_LDADD = $(top_builddir)/libnfft3.la @fftw3_LDFLAGS@ @fftw3_LIBS@

EXTRA_DIST = levelplots.m  show_results.m README
application_DATA = $(EXTRA_DIST)
