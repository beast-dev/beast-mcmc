#!/bin/bash
# exactly the commit the workflow resolved (BEAGLE_SHA), which the cache keys name
git init -q ${BEAGLE_DIR}
cd ${BEAGLE_DIR}
git fetch -q --depth 1 https://github.com/beagle-dev/beagle-lib.git "${BEAGLE_SHA:?BEAGLE_SHA not set}"
git checkout -q FETCH_HEAD
mkdir build
cd build
echo $PWD
cmake -DBUILD_CUDA=OFF -DBUILD_OPENCL=OFF -DBEAGLE_OPTIMIZE_FOR_NATIVE_ARCH=OFF ..
make DESTDIR=${GITHUB_WORKSPACE}/${BEAGLE_DIR} install
#export LD_LIBRARY_PATH=${BEAGLE_LIB}
