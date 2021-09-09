# $Id: Makefile.am 3086 2009-03-11 07:36:56Z keiner $

AM_CPPFLAGS = -I$(top_srcdir)/include

applicationdir = $(pkgdatadir)/applications/fastsumS2

application_PROGRAMS = fastsumS2

fastsumS2_SOURCES = fastsumS2.c
fastsumS2_LDADD = $(top_builddir)/libnfft3.la

EXTRA_DIST = fastsumS2.m writeTestcase.m readTestcase.m README
application_DATA = $(EXTRA_DIST)