#!/bin/bash
# Runs the big-3 comparison: each method with each seed, one after another (so that timings are not
# disturbed by other runs). Each run goes in RUN_DIR/<method>-<seed>/, with BEAST's output in out.txt and
# its wall time (seconds) in time.txt; runs that have finished (time.txt exists) are skipped, so the
# script can be restarted.
#
# usage: run_benchmark.sh XML_DIR RUN_DIR "SEEDS" "METHODS"
#   e.g. run_benchmark.sh xmls runs "11 22" "mat-local mat-path baseline"
#
# BEAST_CP: the classpath of the BEAST build (default: build/dist of this repository plus lib/*)
# JAVA_OPTS: extra JVM options (default: -Xmx8g)
set -u
XML_DIR=$(cd "$1" && pwd)
RUN_DIR=$2
SEEDS=$3
METHODS=$4
REPO=$(cd "$(dirname "$0")/../../.." && pwd)
BEAST_CP=${BEAST_CP:-"$REPO/build/dist/beast.jar:$REPO/lib/*"}
JAVA_OPTS=${JAVA_OPTS:-"-Xmx8g"}
mkdir -p "$RUN_DIR"
RUN_DIR=$(cd "$RUN_DIR" && pwd)

for seed in $SEEDS; do
    for method in $METHODS; do
        dir="$RUN_DIR/$method-$seed"
        if [ -f "$dir/time.txt" ]; then
            echo "$method seed $seed: done already"
            continue
        fi
        mkdir -p "$dir"
        cp "$XML_DIR/$method.xml" "$dir/"
        echo "$method seed $seed: started $(date)"
        start=$(date +%s)
        (cd "$dir" && java $JAVA_OPTS -Dparsers=development -Djava.library.path=/opt/homebrew/lib \
            -cp "$BEAST_CP" dr.app.beast.BeastMain -seed "$seed" -overwrite "$method.xml" > out.txt 2>&1)
        status=$?
        end=$(date +%s)
        if [ $status -ne 0 ]; then
            echo "$method seed $seed: FAILED (status $status), see $dir/out.txt"
            continue
        fi
        echo $((end - start)) > "$dir/time.txt"
        echo "$method seed $seed: finished in $((end - start)) s"
    done
done
