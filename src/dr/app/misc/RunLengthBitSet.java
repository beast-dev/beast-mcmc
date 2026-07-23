/*
 * RunLengthBitSet.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 *  BEAST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 *
 */

package dr.app.misc;

import java.util.Arrays;
import java.util.BitSet;
import java.util.HashSet;
import java.util.Set;

public class RunLengthBitSet {

    private static final int BLOCK_SIZE = 6;

    private int[] runs;
    private int runCount;
    private final int length;

    public RunLengthBitSet(int length) {
        this.length = length;
        runs = new int[BLOCK_SIZE];
        runCount = 0;
    }

    public int getLength() {
        return length;
    }

    /**
     * set a bit (can only set bits not already set)
     * @param index which bit
     */
    public void set(int index) {

        if (index == 0) {
            // easy edge case
            if (runs[0] > 0) {
                runs[0] -= 1;
            }
            runs[1] += 1;
            return;
        }
        // find which run of zeros the bit goes in
        int bitPos = 0;
        int runLen = 0;
        int i = 0; // position in the runs array
        while (index >= bitPos) {
            if (i < runCount) {
                runLen = runs[i];
                bitPos += runLen;
            } else {
                // if we are at the end of the runs array
                // then the remaining bits are zero so
                // skip to the end
                runLen = length - bitPos;
                bitPos = length;
                break;
            }
            if (index >= bitPos) {
                i++; // move to the next run of 1s
                bitPos += runs[i];
                i++; // and move on to the next run of zeros
            }
        }

        if (i == runCount) {
            // the bit being set is at the end
            if (index - (bitPos - runLen)  > 0) {
                insertion(i);
                runs[i] = index - (bitPos - runLen);
                runs[i + 1] = 1;
            } else {
                runs[i - 1] += 1;
            }
        } else {
            if (runLen == 1) {
                // the bit being set is the run of 1 zero
                // merge runs of 1 either side
                runs[i - 1] += runs[i+1] + 1;
                deletion(i);
            } else if (runLen - (bitPos - index) == 0) {
                // we are at the start of the run of zeros
                // so add the bit to the previous run of ones
                runs[i] -= 1;
                runs[i - 1] += 1;
            } else if ((bitPos - index) == 1) {
                // we are at the end of the run of zeros
                // so add the bit to the next run of ones
                runs[i] -= 1;
                runs[i + 1] += 1;
            } else {
                // we are in the middle of a run of zeros
                // so insert a new run
                insertion(i);
                runs[i] -= runLen - (bitPos - index) - 1;
                runs[i + 1] = 1;
                runs[i + 2] += (bitPos - index);
            }
        }
    }

    /**
     * Insert two spaces in runs array, increasing capacity if needed.
     * @param pos position to insert spaces
     */
    private void insertion(int pos) {
        runCount += 2;
        if (runCount > runs.length) {
            // expand array
            int[] newRuns = new int[runs.length + BLOCK_SIZE];
            // Copy everything up to pos...
            // If pos is at the start then this doesn't copy anything
            System.arraycopy(runs, 0, newRuns, 0, pos);
            // ...and then copy everything after pos, shifting up two
            // If pos is at the end then this doesn't copy anything
            System.arraycopy(runs, pos, newRuns, pos + 2, runCount - (pos + 2));
            runs = newRuns;
            assert runCount <= runs.length;
        } else {
            // Just copy everything after pos, shifting up two.
            // If pos is at the end then this doesn't copy anything.
            System.arraycopy(runs, pos, runs, pos + 2, runCount - (pos + 2));
        }
    }

    private void deletion(int pos) {
        System.arraycopy(runs, pos + 2, runs, pos, runCount - (pos + 2));
        runCount -= 2;
        runs[runCount] = 0;
        runs[runCount + 1] = 0; // zero out the unused elements so the arrays can be compared
    }

    public boolean isSet(int index) {
        assert index < length;

        int i = 0;

        int bitPos = 0;
        while (i < runCount) {
            bitPos += runs[i];
            if (index < bitPos) {
                // is a zero
                return false;
            }
            i++;
            bitPos += runs[i];
            if (index < bitPos) {
                // is a one
                return true;
            }
            i++;
        }
        // is a zero
        return false;
    }


    public void union(RunLengthBitSet set) {
        throw new UnsupportedOperationException("not implemented yet");
    }

    @Override
    public boolean equals(Object obj) {
        return Arrays.equals(((RunLengthBitSet)obj).runs, runs);
    }

    @Override
    public int hashCode() {
        return Arrays.hashCode(runs);
    }

    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder();
        int k = 0;
        int i = 0;
        while (i < runCount) {
            for (int j = 0; j < runs[i]; j++) {
                sb.append("0");
            }
            k += runs[i];
            i++;
            for (int j = 0; j < runs[i]; j++) {
                sb.append("1");
            }
            k += runs[i];
            i++;
        }
        for (int j = 0; j < length - k; j++) {
            sb.append("0");
        }

        return sb.toString();
    }

    private static String test(RunLengthBitSet bits) {
        StringBuilder sb = new StringBuilder();
        sb.append(bits);
        sb.append(" - ");

        StringBuilder sb1 = new StringBuilder();
        for (int i = 0; i < bits.runCount; i++) {
            sb1.append(bits.runs[i]);
            sb1.append(" ");
        }
        sb.append(sb1);

        for (int i = 0; i < 15 - sb1.length(); i++) {
            sb.append(" ");
        }

        sb.append(" = ");

        StringBuilder sb2 = new StringBuilder();
        for (int i = 0; i < bits.getLength(); i++) {
            sb2.append(bits.isSet(i) ? "1" : "0");
        }
        sb.append(sb2);

        assert sb2.toString().contentEquals(bits.toString());

        sb.append(" runs count: ");
        sb.append(bits.runCount);
        sb.append(", runs array length: ");
        sb.append(bits.runs.length);
        sb.append(", length: ");
        sb.append(bits.length);

        return sb.toString();
    }

    public static void main(String[] args) {
        int replicates = 500000;
        int length = 1024;
        int bitCount = 256;

        System.out.println("RunLengthBitSet");
        Set<RunLengthBitSet> rlbss = new HashSet<>();
        long startTime = System.currentTimeMillis();
        for (int i = 0; i < replicates; i++) {
            RunLengthBitSet runLengthBitSet = new RunLengthBitSet(length);
            for (int j = 0; j < bitCount; j++) {
                runLengthBitSet.set((int)Math.floor(Math.random() * length));
            }
            rlbss.add(runLengthBitSet);
        }
        System.out.println("Size: " + rlbss.size());

        double elapsed = (double)(System.currentTimeMillis() - startTime) / 1000;
        System.out.println("Time: " + elapsed);

        System.out.println();

        System.out.println("BitSet");
        Set<BitSet> bss = new HashSet<>();
        startTime = System.currentTimeMillis();
        for (int i = 0; i < replicates; i++) {
            BitSet bitSet = new BitSet(length);
            for (int j = 0; j < bitCount; j++) {
                bitSet.set((int)Math.floor(Math.random() * length));
            }
            bss.add(bitSet);
        }
        System.out.println("Size: " + bss.size());
        elapsed = (double)(System.currentTimeMillis() - startTime) / 1000;
        System.out.println("Time: " + elapsed);
    }
}
