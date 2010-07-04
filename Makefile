# these values filled in by    yorick -batch make.i
Y_MAKEDIR=/usr/lib/yorick/2.1
Y_EXE=/usr/lib/yorick/2.1/bin/yorick
Y_EXE_PKGS=
Y_EXE_HOME=/usr/lib/yorick/2.1
Y_EXE_SITE=/usr/share/yorick/2.1

# ----------------------------------------------------- optimization flags

# options for make command line, e.g.-   make COPT=-g TGT=exe
COPT=$(COPT_DEFAULT)
TGT=$(DEFAULT_TGT)

# ------------------------------------------------ macros for this package

PKG_NAME=yao
PKG_I=yao_fast.i yao_utils.i

OBJS=aoSimulUtils.o utils.o yao_fast.o

# change to give the executable a name other than yorick
PKG_EXENAME=yorick

# PKG_DEPLIBS=-Lsomedir -lsomelib   for dependencies of this package
PKG_DEPLIBS=-lfftw3f
# set compiler (or rarely loader) flags specific to this package
PKG_CFLAGS=
PKG_LDFLAGS=

# list of additional package names you want in PKG_EXENAME
# (typically Y_EXE_PKGS should be first here)
EXTRA_PKGS=$(Y_EXE_PKGS)

# list of additional files for clean
PKG_CLEAN=

# autoload file for this package, if any
PKG_I_START=
# non-pkg.i include files for this package, if any
PKG_I_EXTRA=yao.i aoutil.i yaokl.i newfits.i yao_util.i turbulence.i yao_gui.i yaopy.i yao_wfs.i yao_structures.i yao_dm.i yao_svipc.i yao_setnsync.i yaodh.i

# -------------------------------- standard targets and rules (in Makepkg)

# set macros Makepkg uses in target and dependency names
# DLL_TARGETS, LIB_TARGETS, EXE_TARGETS
# are any additional targets (defined below) prerequisite to
# the plugin library, archive library, and executable, respectively
PKG_I_DEPS=$(PKG_I)
Y_DISTMAKE=distmake

include $(Y_MAKEDIR)/Make.cfg
include $(Y_MAKEDIR)/Makepkg
include $(Y_MAKEDIR)/Make$(TGT)

# override macros Makepkg sets for rules and other macros
# Y_HOME and Y_SITE in Make.cfg may not be correct (e.g.- relocatable)
Y_HOME=$(Y_EXE_HOME)
Y_SITE=$(Y_EXE_SITE)

# reduce chance of yorick-1.5 corrupting this Makefile
MAKE_TEMPLATE = protect-against-1.5

# ------------------------------------- targets and rules for this package

# simple example:
#myfunc.o: myapi.h
# more complex example (also consider using PKG_CFLAGS above):
#myfunc.o: myapi.h myfunc.c
#	$(CC) $(CPPFLAGS) $(CFLAGS) -DMY_SWITCH -o $@ -c myfunc.c

clean::
	-rm -rf binaries

install::
	mkdir -p $(DEST_Y_SITE)/python
	mkdir -p $(DEST_Y_SITE)/glade
	mkdir -p $(DEST_Y_SITE)/g
	mkdir -p $(DEST_Y_SITE)/gist
	mkdir -p $(DEST_Y_SITE)/share/yao
	mkdir -p $(DEST_Y_BINDIR)
#	mkdir -p /usr/local/man/man1
	cp -p yao.py $(DEST_Y_SITE)/python/
	cp -p yao.glade $(DEST_Y_SITE)/glade/
	cp -p aosimul3.gs $(DEST_Y_SITE)/g/
	cp -p letter.gs $(DEST_Y_SITE)/g/
	cp -p aosimul3.gs $(DEST_Y_SITE)/gist/
	cp -p letter.gs $(DEST_Y_SITE)/gist/
	cp -p yao $(DEST_Y_BINDIR)/
	cp -pr examples $(DEST_Y_SITE)/share/yao/
	cp -pr doc $(DEST_Y_SITE)/share/yao/
	-rm -rf $(DEST_Y_SITE)/share/yao/examples/CVS
	-rm -rf $(DEST_Y_SITE)/share/yao/doc/CVS
#	gzip -9 doc/yao.1
#	cp -p doc/yao.1 /usr/local/man/man1/

uninstall::
	-rm $(DEST_Y_BINDIR)/yao
	-rm $(DEST_Y_SITE)/g/aosimul3.gs
	-rm $(DEST_Y_SITE)/g/letter.gs
	-rm $(DEST_Y_SITE)/gist/aosimul3.gs
	-rm $(DEST_Y_SITE)/gist/letter.gs
	-rm $(DEST_Y_SITE)/python/yao.py
	-rm $(DEST_Y_SITE)/glade/yao.glade
	-rm -rf $(DEST_Y_SITE)/share/yao/
#	-rm -rf /usr/local/man/man1/yao.1

# -------------------------------------------------------- end of Makefile


# for the binary package production (add full path to lib*.a below):
PKG_DEPLIBS_STATIC=-lm /usr/lib/libfftw3f.a
PKG_ARCH = $(OSTYPE)-$(MACHTYPE)
# The above usually don t work. Edit manually and change the PKG_ARCH below:
# PKG_ARCH = linux-x86
PKG_VERSION = $(shell (awk '{if ($$1=="Version:") print $$2}' $(PKG_NAME).info))
# .info might not exist, in which case he line above will exit in error.

# packages or devel_pkgs:
PKG_DEST_URL = packages

package:
	$(MAKE)
	$(LD_DLL) -o $(PKG_NAME).so $(OBJS) ywrap.o $(PKG_DEPLIBS_STATIC) $(DLL_DEF)
	mkdir -p binaries/$(PKG_NAME)/dist/y_home/lib
	mkdir -p binaries/$(PKG_NAME)/dist/y_home/bin
	mkdir -p binaries/$(PKG_NAME)/dist/y_home/i-start
	mkdir -p binaries/$(PKG_NAME)/dist/y_site/i0
	mkdir -p binaries/$(PKG_NAME)/dist/y_site/g
	mkdir -p binaries/$(PKG_NAME)/dist/y_site/gist
	mkdir -p binaries/$(PKG_NAME)/dist/y_site/python
	mkdir -p binaries/$(PKG_NAME)/dist/y_site/glade
	mkdir -p binaries/$(PKG_NAME)/dist/y_site/share/yao/examples
	mkdir -p binaries/$(PKG_NAME)/dist/y_site/share/yao/doc
	cp -p *.i binaries/$(PKG_NAME)/dist/y_site/i0/
	rm binaries/$(PKG_NAME)/dist/y_site/i0/check.i
	if test -n "$(PKG_I_START)"; then rm binaries/$(PKG_NAME)/dist/y_site/i0/$(PKG_I_START); fi
	cp -p $(PKG_NAME).so binaries/$(PKG_NAME)/dist/y_home/lib/
	if test -f "check.i"; then cp -p check.i binaries/$(PKG_NAME)/.; fi
	if test -n "$(PKG_I_START)"; then cp -p $(PKG_I_START) \
	  binaries/$(PKG_NAME)/dist/y_home/i-start/; fi
	cat $(PKG_NAME).info | sed -e 's/OS:/OS: $(PKG_ARCH)/' > tmp.info
	mv tmp.info binaries/$(PKG_NAME)/$(PKG_NAME).info
	cp -p *.i binaries/$(PKG_NAME)/dist/y_site/i0/.
	cp -p README binaries/$(PKG_NAME)/dist/y_site/share/yao/.
	cp -p INSTALL binaries/$(PKG_NAME)/dist/y_site/share/yao/.
	cp -p LICENSE binaries/$(PKG_NAME)/dist/y_site/share/yao/.
	cp -p yao binaries/$(PKG_NAME)/dist/y_home/bin/.
	cp -p yao.py binaries/$(PKG_NAME)/dist/y_site/python/.
	cp -p yao.glade binaries/$(PKG_NAME)/dist/y_site/glade/.
	cp -p *.gs binaries/$(PKG_NAME)/dist/y_site/g/.
	cp -p *.gs binaries/$(PKG_NAME)/dist/y_site/gist/.
	-cp -p examples/* binaries/$(PKG_NAME)/dist/y_site/share/yao/examples/.
	-cp -p doc/* binaries/$(PKG_NAME)/dist/y_site/share/yao/doc/.
	cd binaries; tar zcvf $(PKG_NAME)-$(PKG_VERSION)-$(PKG_ARCH).tgz $(PKG_NAME)

distbin: package
	if test -f "binaries/$(PKG_NAME)-$(PKG_VERSION)-$(PKG_ARCH).tgz" ; then \
	  ncftpput -f $(HOME)/.ncftp/maumae www/yorick/$(PKG_DEST_URL)/$(PKG_ARCH)/tarballs/ \
	  binaries/$(PKG_NAME)-$(PKG_VERSION)-$(PKG_ARCH).tgz; fi
	if test -f "binaries/$(PKG_NAME)/$(PKG_NAME).info" ; then \
	  ncftpput -f $(HOME)/.ncftp/maumae www/yorick/$(PKG_DEST_URL)/$(PKG_ARCH)/info/ \
	  binaries/$(PKG_NAME)/$(PKG_NAME).info; fi

distsrc:
	make clean; rm -rf binaries
	cd ..; tar --exclude binaries --exclude CVS --exclude .svn --exclude *.spec -zcvf \
	   $(PKG_NAME)-$(PKG_VERSION)-src.tgz yorick-$(PKG_NAME)-$(PKG_VERSION);\
	ncftpput -f $(HOME)/.ncftp/maumae www/yorick/$(PKG_DEST_URL)/src/ \
	   $(PKG_NAME)-$(PKG_VERSION)-src.tgz
	ncftpput -f $(HOME)/.ncftp/maumae www/yorick/contrib/ \
	   ../$(PKG_NAME)-$(PKG_VERSION)-src.tgz


# -------------------------------------------------------- end of Makefile
