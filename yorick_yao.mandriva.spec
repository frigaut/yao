%define name yorick-yao
%define version 4.8.2
%define release gemini2010sep16

Summary: yorick adaptive optics simulation package
Name: %{name}
Version: %{version}
Release: %{release}
Source0: %{name}-%{version}.tar.bz2
License: BSD
Group: Applications/Engineering
Packager: Francois Rigaut <frigaut@gemini.edu>
Url: http://www.maumae.net/yorick/doc/plugins.php
BuildRoot: %{_tmppath}/%{name}-%{version}-%{release}-buildroot
Requires: yorick >= 2.1 yorick-imutil >= 0.5 yorick-yutils >= 1.0 libfftw3


%description
Multi-usage adaptive optics Monte-Carlo simulation package.
Can simulate versatile configs (classical, LGS, GLAO, MCAO)
with Shack-Hartmann or curvature sensors. Fast.

From version 4.0, yao has now a gtk interface.

%prep
%setup -q

%build
yorick -batch make.i
make
if [ -f check.i ] ; then
   mv check.i %{name}_check.i
fi;

%install
rm -rf $RPM_BUILD_ROOT
mkdir -p $RPM_BUILD_ROOT/usr/lib/yorick/lib
mkdir -p $RPM_BUILD_ROOT/usr/lib/yorick/i0
mkdir -p $RPM_BUILD_ROOT/usr/lib/yorick/i
mkdir -p $RPM_BUILD_ROOT/usr/lib/yorick/g
mkdir -p $RPM_BUILD_ROOT/usr/lib/yorick/glade
mkdir -p $RPM_BUILD_ROOT/usr/lib/yorick/python
mkdir -p $RPM_BUILD_ROOT/usr/share/doc/yorick-yao
mkdir -p $RPM_BUILD_ROOT/usr/share/man/man1
mkdir -p $RPM_BUILD_ROOT/usr/bin
mkdir -p $RPM_BUILD_ROOT/usr/lib/yorick/packages/installed

install -m 755 yao.so $RPM_BUILD_ROOT/usr/lib/yorick/lib
install -m 755 yao $RPM_BUILD_ROOT/usr/bin
install -m 644 yao_fast.i $RPM_BUILD_ROOT/usr/lib/yorick/i0
install -m 644 yao_utils.i $RPM_BUILD_ROOT/usr/lib/yorick/i0
install -m 644 *.i $RPM_BUILD_ROOT/usr/lib/yorick/i
install -m 644 *.gs $RPM_BUILD_ROOT/usr/lib/yorick/g
install -m 644 yao.glade $RPM_BUILD_ROOT/usr/lib/yorick/glade
install -m 755 yao.py $RPM_BUILD_ROOT/usr/lib/yorick/python
install -m 644 LICENSE $RPM_BUILD_ROOT/usr/share/doc/yorick-yao
install -m 644 README $RPM_BUILD_ROOT/usr/share/doc/yorick-yao
install -dm 755 examples $RPM_BUILD_ROOT/usr/share/doc/yorick-yao/examples
install -m 644 examples/*.par $RPM_BUILD_ROOT/usr/share/doc/yorick-yao/examples
install -m 644 examples/*.i $RPM_BUILD_ROOT/usr/share/doc/yorick-yao/examples
install -m 755 examples/testclean $RPM_BUILD_ROOT/usr/share/doc/yorick-yao/examples
install -m 644 doc/yao.1.gz $RPM_BUILD_ROOT/usr/share/man/man1
install -m 644 yao.info $RPM_BUILD_ROOT/usr/lib/yorick/packages/installed

rm $RPM_BUILD_ROOT/usr/lib/yorick/i/yao_fast.i
rm $RPM_BUILD_ROOT/usr/lib/yorick/i/yao_utils.i

%clean
rm -rf $RPM_BUILD_ROOT

%files
%defattr(-,root,root)
/usr/lib/yorick/lib/yao.so
/usr/bin/yao
/usr/lib/yorick/i0/*.i
/usr/lib/yorick/i/*.i
/usr/lib/yorick/g/*.gs
/usr/lib/yorick/glade/yao.glade
/usr/lib/yorick/python/yao.py
/usr/share/doc/yorick-yao
/usr/share/man/man1
/usr/lib/yorick/packages/installed/*

%changelog
* Tue Jan 09 2008 <frigaut@users.sourceforge.net>
- included the info file for compat with pkg_mngr

* Mon Dec 31 2007 <frigaut@users.sourceforge.net>
- new distro directory structure
- updated to cvs
