#!/bin/bash
git clone -b ${BEAGLE_BRANCH} --depth 1 https://github.com/beagle-dev/beagle-lib.git ${BEAGLE_DIR}
cd ${BEAGLE_DIR}
mkdir build
cd build
echo $PWD
cmake -DBUILD_CUDA=OFF -DBUILD_OPENCL=OFF -DBEAGLE_OPTIMIZE_FOR_NATIVE_ARCH=OFF ..
make DESTDIR=${GITHUB_WORKSPACE}/${BEAGLE_DIR} install
#export LD_LIBRARY_PATH=${BEAGLE_LIB}
