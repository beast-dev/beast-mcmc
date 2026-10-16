#!/bin/bash
git clone -b ${BEAGLE_BRANCH} --depth 1 https://github.com/beagle-dev/beagle-lib.git ${BEAGLE_DIR}
cd ${BEAGLE_DIR}
# v4_release uses set(), not option(), so the -D flag below is ignored; force it off.
sed -i 's/set(BEAGLE_OPTIMIZE_FOR_NATIVE_ARCH true)/set(BEAGLE_OPTIMIZE_FOR_NATIVE_ARCH false)/' CMakeLists.txt
mkdir build
cd build
cmake -DBUILD_CUDA=OFF -DBUILD_OPENCL=OFF -DBEAGLE_OPTIMIZE_FOR_NATIVE_ARCH=OFF ..
make DESTDIR=${GITHUB_WORKSPACE}/${BEAGLE_DIR} install
#export LD_LIBRARY_PATH=${BEAGLE_LIB}
