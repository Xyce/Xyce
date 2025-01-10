SUBDIRS = XyceCInterface

bin_SCRIPTS = $(MAYBE_SCRIPTS)
EXTRA_SCRIPTS=buildxyceplugin
CLEANFILES = $(bin_SCRIPTS)
transform = @program_transform_name@

do_subst = sed -e 's|[@]datadir[@]|$(datadir)|g' \
	-e 's|[@]bindir[@]|$(bindir)|g' \
	-e 's|[@]datadir[@]|$(datadir)|g' \
	-e 's|[@]includedir[@]|$(includedir)|g' \
	-e 's|[@]libexecdir[@]|$(libexecdir)|g' \
	-e 's|[@]libdir[@]|$(libdir)|g' \
	-e 's|[@]INSTALL[@]|$(INSTALL)|g' \
	-e 's|[@]RAW_CPPFLAGS[@]|$(RAW_CPPFLAGS)|g' \
	-e 's|[@]CXXFLAGS[@]|$(CXXFLAGS)|g' \
	-e 's|[@]LDFLAGS[@]|$(LDFLAGS)|g' \
	-e 's|[@]CXX[@]|$(CXX)|g' \
	-e 's|[@]LIBTOOL[@]|'`echo libtool | sed '$(transform)'`'|g'

buildxyceplugin: buildxyceplugin.in Makefile
	$(do_subst) < $(srcdir)/buildxyceplugin.in  > buildxyceplugin
	chmod +x buildxyceplugin
