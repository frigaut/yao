%define name yorick_yao
%define version 4.8.2
%define release 01gemini
%define yorick_version 2.2

Summary: yorick adaptive optics simulation package
Name: %{name}
Version: %{version}
Release: %{release}
Source0: %{name}-%{version}.tar.bz2
License: BSD
Group: Development/Languages Applications/Engineering
Packager: Francois Rigaut <frigaut@gemini.edu>
Url: http://www.maumae.net/yorick/doc/plugins.php
BuildRoot: %{_tmppath}/%{name}-%{version}-%{release}-buildroot
Requires: yorick = %{yorick_version} yorick_imutil >= 0.5 yorick_utils >= 1.0 libfftw3


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
mkdir -p $RPM_BUILD_ROOT/usr/lib/yorick/%{yorick_version}/lib
mkdir -p $RPM_BUILD_ROOT/usr/lib/yorick/%{yorick_version}/bin
mkdir -p $RPM_BUILD_ROOT/usr/lib/yorick/%{yorick_version}/i0
mkdir -p $RPM_BUILD_ROOT/usr/lib/yorick/%{yorick_version}/g
mkdir -p $RPM_BUILD_ROOT/usr/lib/yorick/%{yorick_version}/glade
mkdir -p $RPM_BUILD_ROOT/usr/lib/yorick/%{yorick_version}/python
mkdir -p $RPM_BUILD_ROOT/usr/share/doc
mkdir -p $RPM_BUILD_ROOT/usr/bin

install -dm 755 examples $RPM_BUILD_ROOT/usr/share/yao/examples
install -m 644 examples/*.par $RPM_BUILD_ROOT/usr/share/yao/examples
install -m 644 examples/*.i $RPM_BUILD_ROOT/usr/share/yao/examples
install -m 755 examples/testclean $RPM_BUILD_ROOT/usr/share/yao/examples
install -m 755 yao.so $RPM_BUILD_ROOT/usr/lib/yorick/%{yorick_version}/lib
install -m 644 *.i $RPM_BUILD_ROOT/usr/lib/yorick/%{yorick_version}/i0
install -m 644 *.gs $RPM_BUILD_ROOT/usr/lib/yorick/%{yorick_version}/g
install -m 644 yao.glade $RPM_BUILD_ROOT/usr/lib/yorick/%{yorick_version}/glade
install -m 755 yao.py $RPM_BUILD_ROOT/usr/lib/yorick/%{yorick_version}/python
install -m 644 README $RPM_BUILD_ROOT/usr/share/yao
install -m 755 yaogtk $RPM_BUILD_ROOT/usr/bin

%clean
rm -rf $RPM_BUILD_ROOT

%files
%defattr(-,root,root)
/usr/lib/yorick/%{yorick_version}/lib/yao.so
/usr/bin/yaogtk
/usr/lib/yorick/%{yorick_version}/i0/*.i
/usr/lib/yorick/%{yorick_version}/g/*.gs
/usr/lib/yorick/%{yorick_version}/glade/yao.glade
/usr/lib/yorick/%{yorick_version}/python/yao.py
/usr/share/yao


%changelog
