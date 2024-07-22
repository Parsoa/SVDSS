#!/bin/bash

set -ex

# Activate Holy Build Box environment.
source /hbb_exe/activate

# Enter source
cd io

# Install deps
# This doesn't work
# yum install -y zlib-devel bzip2-devel xz-devel
# since I want static libs. And -devel packages do not provide static .a.
# This is a workaround to get them from the rpm packages

# yumdownloader --source bzip2-devel # this doesn't work anymore (due to centos7 eol?)
# as long as the package is hosted here, I can do this
curl https://vault.centos.org/7.9.2009/os/Source/SPackages/bzip2-1.0.6-13.el7.src.rpm > bzip2-1.0.6-13.el7.src.rpm
rpm2cpio bzip2-1.0.6-13.el7.src.rpm | cpio -idv
tar xvfz bzip2-1.0.6.tar.gz
cd bzip2-1.0.6
make libbz2.a
ls libbz2.a
cd ..

# yumdownloader --source xz-devel
curl https://vault.centos.org/7.9.2009/os/Source/SPackages/xz-5.2.2-1.el7.src.rpm > xz-5.2.2-1.el7.src.rpm
rpm2cpio xz-5.2.2-1.el7.src.rpm | cpio -idv
tar xvfz xz-5.2.2.tar.gz
cd xz-5.2.2
./configure --enable-static
make
ls src/liblzma/.libs/liblzma.a
cd ..

# Compile
rm -rf holy-build
mkdir -p holy-build
cd holy-build
cmake -DCMAKE_BUILD_TYPE=Release -DHOLYBUILD=ON ..
make

# Check and copy result to host
# hardening-check -b ../SVDSS
libcheck ../SVDSS
mv ../SVDSS /io/holy/SVDSS_linux_x86-64

# clean
rm -rf /io/bzip2* /io/xz*
