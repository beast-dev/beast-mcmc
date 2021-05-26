#!/bin/bash

passed=true
failedFiles=()

for file in ci/TestXML/*\.xml
do
  if java -Djava.library.path=${BEAGLE_LIB} -jar build/dist/beast.jar -fail_threads -seed 666 -overwrite $file; then
    echo $file passed
  else
    echo $file failed
    passed=false
    failedFiles+=($file)
  fi
done

echo
echo

if $passed
then
  echo "All test xml passed."
else
  echo "The following test xml failed:"
  for file in "${failedFiles[@]}"
  do
    echo "$file failed"
  done
  exit -1
fi