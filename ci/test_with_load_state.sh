#!/bin/bash
cd ci

BEAGLE_FILES=${GITHUB_WORKSPACE}/${BEAGLE_LIB}
echo $BEAGLE_FILES
ls $BEAGLE_FILES

java -jar -Djava.library.path=${BEAGLE_FILES} build/dist/beast.jar -beagle_info


for file in TestXMLwithLoadState/*\.xml
do
  filename=$(basename $file .xml)
  if java -Djava.library.path=$BEAGLE_FILES -jar ../build/dist/beast.jar -fail_threads -seed 666 -load_state $filename.chkpt -overwrite $file; then
    echo $file passed
  else
    echo $file failed; exit -1
  fi
done
