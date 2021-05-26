#!/bin/bash
git clone -b ${BEAGLE_BRANCH} --depth 1 https://github.com/beagle-dev/beagle-lib.git ${BEAGLE_DIR}
cd ${BEAGLE_DIR}
mkdir build
cd build
echo $PWD
cmake -DBUILD_CUDA=OFF -DBUILD_OPENCL=OFF ..
make DESTDIR=${BEAGLE_DIR} install
export LD_LIBRARY_PATH=${BEAGLE_LIB}
cd ${BEAGLE_LIB}
ls