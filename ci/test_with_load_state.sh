#!/bin/bash

java -jar -Djava.library.path=${BEAGLE_LIB} build/dist/beast.jar -beagle_info


for file in ci/TestXMLwithLoadState/*\.xml
do
  checkpoint=ci/TestXMLwithLoadState/$(basename $file .xml).chkpt
  if java -Djava.library.path=${BEAGLE_LIB} -jar build/dist/beast.jar -fail_threads -seed 666 -load_state $checkpoint -overwrite $file; then
    echo $file passed
  else
    echo $file failed; exit -1
  fi
done
