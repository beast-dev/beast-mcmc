#!/bin/bash
git clone -b ${BEAGLE_BRANCH} --depth 1 https://github.com/beagle-dev/beagle-lib.git ${BEAGLE_DIR}
cd ${BEAGLE_DIR}
# v4_release ignores -DBEAGLE_OPTIMIZE_FOR_NATIVE_ARCH (plain set(), not option()),
# and its SSE plugin needs SSE4.1, so swap -march=native for a portable baseline
# instead of removing it.
sed -i 's/-march=native/-march=x86-64-v2/g' CMakeLists.txt
mkdir build
cd build
cmake -DBUILD_CUDA=OFF -DBUILD_OPENCL=OFF -DBEAGLE_OPTIMIZE_FOR_NATIVE_ARCH=OFF ..
make DESTDIR=${GITHUB_WORKSPACE}/${BEAGLE_DIR} install
#export LD_LIBRARY_PATH=${BEAGLE_LIB}
