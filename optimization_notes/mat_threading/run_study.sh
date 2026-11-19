#!/bin/bash
# Runs the whole profiling study: every XML with the profiles on, a subset with them off (to measure their
# overhead), and a subset with a flight recording. Runs are sequential; finished ones are skipped.
#
# usage: run_study.sh XML_DIR RUNS_DIR
set -u
XML_DIR=$1
RUNS_DIR=$2
HERE=$(cd "$(dirname "$0")" && pwd)
ALL="small-2-local small-2-path big-3-local big-3-path big-3-gamma-local big-2-local big-2-path wide-1-local wide-1-path sim-20k-r005 sim-20k-r02 sim-100k-r005 sim-100k-r02 sim-300k-r02"
OFF="small-2-local big-3-local big-2-local sim-100k-r02"
JFR_RUNS="big-3-local big-2-local sim-100k-r02 sim-300k-r02"
export JAVA_OPTS=${JAVA_OPTS:-"-Xmx16g"}
PROFILE=1 "$HERE/run_profile.sh" "$XML_DIR" "$RUNS_DIR/profiled" "$ALL"
PROFILE=0 "$HERE/run_profile.sh" "$XML_DIR" "$RUNS_DIR/unprofiled" "$OFF"
PROFILE=1 JFR=1 "$HERE/run_profile.sh" "$XML_DIR" "$RUNS_DIR/jfr" "$JFR_RUNS"
echo "study finished $(date)"
