#!/bin/bash

set -e

cd "$(dirname "$0")/.."

XML="tests/BenchmarkXML/canonical_ou_orthogonal_hmc_gradient_check.xml"
LOG="tests/BenchmarkXML/canonical_ou_orthogonal_hmc_gradient_check.log"

java -Djava.library.path="${BEAGLE_LIB}" \
  -cp "build:lib/*" \
  dr.app.beast.BeastMain \
  -overwrite \
  -seed 666 \
  "${XML}"

rm -f "${LOG}"
