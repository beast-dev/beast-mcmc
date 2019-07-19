#!/bin/bash
for file in TestXML/*\.xml
do
  if java -jar ../build/dist/beast.jar -overwrite $file; then
    echo $file passed
  else
    echo $file failed; exit -1
  fi
done
