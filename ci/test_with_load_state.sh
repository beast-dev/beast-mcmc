#!/bin/bash
for file in TestXMLwithLoadState/*\.xml
do
  filename=$(basename $file .xml)
  if java -Djava.library.path=${BEAGLE_LIB} -jar build/dist/beast.jar -fail_threads -seed 666 -load_state $filename.chkpt -overwrite $file; then
    echo $file passed
  else
    echo $file failed; exit -1
  fi
done
