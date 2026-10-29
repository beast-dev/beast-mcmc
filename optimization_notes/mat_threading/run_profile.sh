#!/bin/bash
# Runs the profiling study's XMLs one after another with the step and MAT profiles on, each in
# RUN_DIR/<name>/ with BEAST's output in out.txt and the wall time (seconds) in time.txt; finished runs
# (time.txt exists) are skipped. With JFR=1 a flight recording (run.jfr) is made as well.
#
# usage: run_profile.sh XML_DIR RUN_DIR "NAMES" [SEED]
#   e.g. run_profile.sh xmls runs "big-3-local big-3-path" 11
#
# BEAST_CP: the classpath of the BEAST build (default: build/dist/beast.jar of this repository plus lib/*; make it with ant dist)
# JAVA_OPTS: extra JVM options (default: -Xmx8g)
# PROFILE: 1 (default) for -Dmcmc.profile=true -Dmat.profile=true, 0 for a run without the timers
set -u
XML_DIR=$(cd "$1" && pwd)
RUN_DIR=$2
NAMES=$3
SEED=${4:-11}
REPO=$(cd "$(dirname "$0")/../.." && pwd)
BEAST_CP=${BEAST_CP:-"$REPO/build/dist/beast.jar:$REPO/lib/*"}
JAVA_OPTS=${JAVA_OPTS:-"-Xmx8g"}
PROFILE=${PROFILE:-1}
JFR=${JFR:-0}
mkdir -p "$RUN_DIR"
RUN_DIR=$(cd "$RUN_DIR" && pwd)

for name in $NAMES; do
    dir="$RUN_DIR/$name"
    if [ -f "$dir/time.txt" ]; then
        echo "$name: done already"
        continue
    fi
    mkdir -p "$dir"
    cp "$XML_DIR/$name.xml" "$dir/"
    opts="$JAVA_OPTS"
    if [ "$PROFILE" = 1 ]; then
        opts="$opts -Dmcmc.profile=true -Dmat.profile=true"
    fi
    if [ "$JFR" = 1 ]; then
        opts="$opts -XX:StartFlightRecording=filename=$dir/run.jfr,settings=profile"
    fi
    echo "$name: started $(date)"
    start=$(date +%s)
    (cd "$dir" && java $opts -Dparsers=development -Djava.library.path=/opt/homebrew/lib \
        -cp "$BEAST_CP" dr.app.beast.BeastMain -seed "$SEED" -overwrite "$name.xml" > out.txt 2>&1)
    status=$?
    end=$(date +%s)
    if [ $status -ne 0 ]; then
        echo "$name: FAILED (status $status), see $dir/out.txt"
        continue
    fi
    echo $((end - start)) > "$dir/time.txt"
    echo "$name: finished in $((end - start)) s"
done
