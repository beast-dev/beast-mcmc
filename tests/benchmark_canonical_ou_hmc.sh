#!/bin/bash
# Benchmark of the canonical OU HMC on the 81-taxon, 3-trait orthogonal-block example.
#
# The data and model live in one place, tests/TestXML/testCanonicalOUOrthogonalHmcGradientCheck.xml
# (a 3-state smoke test).  This script derives the benchmark run from it: 1000 HMC states,
# gradient checks on, screen log every 100 states.

set -e

cd "$(dirname "$0")/.."

SRC="tests/TestXML/testCanonicalOUOrthogonalHmcGradientCheck.xml"
WORK="$(mktemp -d)"
XML="${WORK}/canonical_ou_orthogonal_hmc_gradient_check.xml"
trap 'rm -rf "${WORK}"' EXIT

sed -e 's/chainLength="3"/chainLength="1000"/' \
    -e 's/gradientCheckCount="3"/gradientCheckCount="1000"/' \
    -e 's/logEvery="1"/logEvery="100"/' \
    "${SRC}" > "${XML}"

java -Djava.library.path="${BEAGLE_LIB}" \
  -cp "build:lib/*" \
  dr.app.beast.BeastMain \
  -overwrite \
  -seed 666 \
  "${XML}"
