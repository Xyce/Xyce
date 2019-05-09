# Process this file with "automake --foreign" to produce a Makefile.in
# NEVER edit the Makefile.in produced by automake or the Makefile produced
# by configure, as it will get overwritten when automake or configure
# are run

ACLOCAL_AMFLAGS=-I config
LIBTOOL_DEPS = @LIBTOOL_DEPS@

# define which subdirectories to do recursive makes in.
SUBDIRS=src $(MAYBE_PLUGIN) $(MAYBE_UTILS)
libexec_SCRIPTS=$(MAYBE_LIBEXEC_SCRIPTS)
EXTRA_SCRIPTS=libtool

if BUILD_XYCE_SHAREABLE

plugin:
	echo 'Making plugin in user_plugin'
	cd user_plugin && $(MAKE) $(AM_MAKEFLAGS) plugin

install-plugin:
	echo 'Making install-plugin in user_plugin'
	cd user_plugin && $(MAKE) $(AM_MAKEFLAGS) install-plugin

endif

dist-hook::
	rm -rf `find $(top_distdir) -type d -name .svn`
	rm -rf `find $(top_distdir) -type d -name .git`
	rm -rf `find $(top_distdir) -type d -name CVS`

libtool: $(LIBTOOL_DEPS)
	$(SHELL) ./config.status libtool
