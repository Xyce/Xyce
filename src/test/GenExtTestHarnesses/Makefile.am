
AM_CPPFLAGS = @Xyce_INCS@

# standalone executable
check_PROGRAMS = testGenCoup TestGH19
testGenCoup_SOURCES = testGenCoup.C
testGenCoup_LDADD = $(top_builddir)/src/libxyce.la
testGenCoup_LDFLAGS = -static $(AM_LDFLAGS)

TestGH19_SOURCES = TestGH19.C
TestGH19_LDADD = $(top_builddir)/src/libxyce.la
TestGH19_LDFLAGS = -static $(AM_LDFLAGS)
