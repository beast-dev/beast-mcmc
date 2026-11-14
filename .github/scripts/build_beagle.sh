#!/bin/bash
set -e
# exactly the commit the workflow resolved (BEAGLE_SHA), which the cache keys name
git init -q ${BEAGLE_DIR}
cd ${BEAGLE_DIR}
git fetch -q --depth 1 https://github.com/beagle-dev/beagle-lib.git "${BEAGLE_SHA:?BEAGLE_SHA not set}"
git checkout -q FETCH_HEAD
# no -march=native, so a cached library runs on every runner: BEAGLE 4.2 takes -DBEAGLE_OPTIMIZE_FOR_NATIVE_ARCH=OFF,
# older releases hard-code set(BEAGLE_OPTIMIZE_FOR_NATIVE_ARCH true)
perl -0pi -e 's/set\(BEAGLE_OPTIMIZE_FOR_NATIVE_ARCH true\)/set(BEAGLE_OPTIMIZE_FOR_NATIVE_ARCH false)/' CMakeLists.txt
grep -q -E 'set\(BEAGLE_OPTIMIZE_FOR_NATIVE_ARCH false\)|option\(BEAGLE_OPTIMIZE_FOR_NATIVE_ARCH' CMakeLists.txt
mkdir build
cd build
echo $PWD
# -msse4.1: BEAGLE before 4.2 needs SSE4.1 (_mm_dp_pd) in its SSE plugin, which -march=native used to provide
# -DBUILD_BIT=ON: the BASTA JNI library (libhmsbeagle-jni-bit) that BEAST's BastaJNIWrapper loads; off by default in
# 4.1 and 4.2, ignored by 4.0.x (which has no BASTA)
cmake -DBUILD_CUDA=OFF -DBUILD_OPENCL=OFF -DBEAGLE_OPTIMIZE_FOR_NATIVE_ARCH=OFF -DBUILD_BIT=ON -DCMAKE_CXX_FLAGS="-msse4.1" ..
make DESTDIR=${GITHUB_WORKSPACE}/${BEAGLE_DIR} install
#export LD_LIBRARY_PATH=${BEAGLE_LIB}
