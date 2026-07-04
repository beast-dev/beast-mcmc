#!/bin/bash
set -e
git clone -b ${BEAGLE_BRANCH} --depth 1 https://github.com/beagle-dev/beagle-lib.git ${BEAGLE_DIR}
cd ${BEAGLE_DIR}
perl -0pi -e 's/set\(BEAGLE_OPTIMIZE_FOR_NATIVE_ARCH true\)/set(BEAGLE_OPTIMIZE_FOR_NATIVE_ARCH false)/' CMakeLists.txt
grep -q 'set(BEAGLE_OPTIMIZE_FOR_NATIVE_ARCH false)' CMakeLists.txt
mkdir build
cd build
echo $PWD
cmake -DBUILD_CUDA=OFF -DBUILD_OPENCL=OFF ..
make DESTDIR=${GITHUB_WORKSPACE}/${BEAGLE_DIR} install
#export LD_LIBRARY_PATH=${BEAGLE_LIB}
