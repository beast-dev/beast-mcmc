#!/usr/bin/env python3
"""
Summarises a Java Flight Recorder file: the top frames of the execution samples (self time) and the
frames most often anywhere on the stack (inclusive time), as counts and percentages.

usage: jfr_top.py RECORDING.jfr [--top 40] [--depth 64]
(needs the `jfr` tool of the JDK on the PATH, or JAVA_HOME set)
"""
import argparse
import collections
import os
import re
import subprocess
import sys


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('recording')
    parser.add_argument('--top', type=int, default=40)
    parser.add_argument('--depth', type=int, default=64)
    parser.add_argument('--prefix', default='dr.', help='only count frames whose class starts with this (inclusive counts)')
    args = parser.parse_args()
    jfr = 'jfr'
    if os.environ.get('JAVA_HOME'):
        jfr = os.path.join(os.environ['JAVA_HOME'], 'bin', 'jfr')
    text = subprocess.run([jfr, 'print', '--events', 'jdk.ExecutionSample', '--stack-depth', str(args.depth),
                           args.recording], capture_output=True, text=True, check=True).stdout
    self_counts = collections.Counter()
    inclusive = collections.Counter()
    samples = 0
    for event in text.split('jdk.ExecutionSample'):
        frames = re.findall(r'^\s+([\w.$<>]+\([^)]*\))(?:\s+line: \d+)?\s*$', event, re.M)
        if not frames:
            continue
        samples += 1
        self_counts[frames[0]] += 1
        for frame in set(frames):
            if frame.startswith(args.prefix):
                inclusive[frame] += 1
    print('%d samples' % samples)
    print('\nSelf time (top frame):')
    for frame, count in self_counts.most_common(args.top):
        print('%6.2f%%  %s' % (100.0 * count / samples, frame))
    print('\nInclusive (anywhere on the stack, %s*):' % args.prefix)
    for frame, count in inclusive.most_common(args.top):
        print('%6.2f%%  %s' % (100.0 * count / samples, frame))


if __name__ == '__main__':
    main()
