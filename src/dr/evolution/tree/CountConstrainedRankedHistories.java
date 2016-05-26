/*
 * CountConstrainedRankedHistories.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
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
 */

package dr.evolution.tree;

import dr.math.Binomial;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.List;

/**
 * @author Alexei Drummond
 */
public class CountConstrainedRankedHistories {


    private static BitSet createClade(String s) {
        BitSet bitSet = new BitSet();
        for (int i = 0; i < s.length(); i++) {
            if (s.charAt(i) == '1') {
                bitSet.set(i);
            }
        }
        return bitSet;
    }

    /**
     * Compute the number of 'free' nodes under each calibration
     *
     * @param constraints
     * @param f
     */
    private static void computeF(BitSet[] constraints, int[] f) {

        for (int i = 0; i < f.length; i++) {

            BitSet set = new BitSet();
            set.or(constraints[i]);

            int subCladeCount = 0;

            for (int j = i + 1; j < f.length; j++) {

                int parent = parent(constraints, j);

                if (parent == i) {
                    subCladeCount += 1;
                    set.andNot(constraints[j]);
                }
            }
            f[i] = set.cardinality() + subCladeCount - 2;
        }
    }

    private static long Nc(int[] f) {

        int C = f.length;
        long Nc = 1;
        for (int i = 0; i < f.length; i++) {
            if (f[i] != 0) Nc *= Binomial.choose(f[i] + C - i - 1, C - i - 1);
        }
        return Nc;
    }

    // returns the number of ranked forests with n tips and k internal nodes

    private static long R(int n, int k) {
        long Rnk = 1;
        for (int i = 0; i < k; i++) {
            Rnk *= Binomial.choose2(n - i);
        }

        if (Rnk == 0) {
            //System.out.println("R(" + n + ", " + k + ") = " + 0);

        }

        return Rnk;
    }

    private static long X(int[][] n, int[][] k, int[] f) {

        int C = f.length;

        long X = 1;

        // for all levels
        for (int i = 0; i < C; i++) {
            // for all forests in level i
            for (int j = 0; j < i + 1; j++) {

//                if (n[i][j] == 1 && k[i][j] == 1) {
//                    System.out.println("n(" + i + ", " + j + ") == 1 && k(" + i + ", " + j + ") == 1");
//                }

                X *= R(n[i][j], k[i][j]);
            }
        }

        // for all levels
        for (int i = 0; i < C; i++) {
            X *= totalOrder(k[i]);
        }
        return X;
    }

    private static int[] nC(BitSet[] clades) {

        int[] nC = new int[clades.length];

        for (int i = 0; i < clades.length; i++) {
            nC[i] = clades[i].cardinality();

            BitSet set = new BitSet();
            set.or(clades[i]);


            for (int j = i + 1; j < clades.length; j++) {
                if (clades[j].intersects(clades[i])) {
                    set.andNot(clades[j]);
                }
            }
            nC[i] = set.cardinality();
        }
        return nC;
    }


    // return the number of total orders of k.length partial orders of lengths k[0], k[1], ..., k[k.length-1]

    private static long totalOrder(int[] k) {

        long T = 1;

        int m = k[0];
        int n = 0;
        for (int i = 0; i < k.length - 1; i++) {
            m += k[i + 1];
            n += k[i];
            T *= Binomial.choose(m, n);
        }
        return T;
    }

    private static int parent(BitSet[] clades, int child) {

        for (int parent = child - 1; parent >= 0; parent--) {
            if (clades[parent].intersects(clades[child])) {
                return parent;
            }
        }
        return 0;
    }

    private static void computeN(int[][] k, int[][] n, int[] nC, BitSet[] clades) {

        int C = nC.length;

        for (int i = 0; i < C; i++) {
            n[C - 1][i] = nC[i];
        }

        for (int level = C - 1; level > 0; level--) {
            for (int clade = 0; clade < level; clade++) {
                n[level - 1][clade] = n[level][clade] - k[level][clade];

            }
        }

        for (int level = C - 1; level > 0; level--) {
            int parent = parent(clades, level);

            //System.out.println("parent = " + parent);

            for (int c = level - 1; c >= parent; c--) {
                n[c][parent] += 1;
            }
        }
    }

    private static boolean kton_next(int[][] v, int[] n, int level) {

        int k = v[level].length - 1;

        if (k < 0) {
            if (level < v.length-1) {

                reset(v, level);

                return kton_next(v, n, level + 1);
                //return true;
            } else return false;
        }

        //System.out.println("level = " + level + ", k = " + k);

        while (v[level][k] == n[level] - 1) {
            k -= 1;
            if (k < 0) {
                // reset and advance 1 at the next level up, or if at top level then end

                if (level < v.length-1) {

                    reset(v, level);
                    return kton_next(v, n, level + 1);
                    //return true;
                } else return false;
            }
        }
        int vk = v[level][k] + 1;

        while (k < v[level].length) {
            v[level][k] = vk;
            k += 1;
        }
        return true;
    }

    private static void reset(int[][] v, int level) {
        for (int l = 0; l < level; l++) {
            for (int i = 0; i < v[l].length; i++) {
                v[l][i] = 0;
            }
        }
    }

    public static void print(int[][] k, int[][] n) {

        System.out.println("k:");
        for (int level = 0; level < k.length; level++) {
            for (int sublevel = 0; sublevel < k[level].length; sublevel++) {

                int subforest = sublevel + level;

                System.out.print(k[level][sublevel] + "\t");

            }
            System.out.println();
        }

        System.out.println("n:");
        for (int level = 0; level < k.length; level++) {
            for (int sublevel = 0; sublevel < k[level].length; sublevel++) {
                System.out.print(n[level][sublevel] + "\t");

            }
            System.out.println();
        }

    }


    public static void main(String[] args) {

        //String[] constraintStrings = {"1111111", "1110000", "0000011"};

        String[] constraintStrings = {
                "11111111111111111111",
                "00000000011110000000",
                "00000000000110000000",
                "00000000011000000000",
                "00000000000000011111",
                "11111111100000000000",
                "11111111000000000000",
                "11111000000000000000"};
//
//        String[] constraintStrings = {
//                "111111111",
//                "111111110",
//                "111110000"};


        // number of constraints, including the "1,1,...,1,1" constraint
        int C = constraintStrings.length;

        System.out.println("n = " + constraintStrings[0].length());
        System.out.println("#calibrations = " + (constraintStrings.length - 1));

        System.out.println("Constraints:");
        for (int i = 0; i < C; i++) {
            System.out.println("  " + constraintStrings[i]);
        }

        BitSet[] constraints = new BitSet[C];

        for (int i = 0; i < C; i++) {
            constraints[i] = createClade(constraintStrings[i]);
        }

        int[] f = new int[C];

        computeF(constraints, f);

//        for (int i = 0; i < f.length; i++) {
//            System.out.println((i + 1) + " : " + f[i]);
//        }

        //System.out.println("Nc = " + Nc(f));

        int[] nC = nC(constraints);
//        for (int i = 0; i < nC.length; i++) {
//            System.out.println("n_{" + (C - 1) + "," + i + "} = " + nC[i]);
//        }


        //System.out.println("T(1,2,1) = " + totalOrder(new int[]{1, 2, 1}));

        int[][] k = new int[C][];
        int[][] n = new int[C][];
        for (int i = 0; i < k.length; i++) {
            k[i] = new int[i + 1];
            n[i] = new int[i + 1];
        }

        computeN(k, n, nC, constraints);

//        for (int i = 0; i < n.length; i++) {
//            for (int j = 0; j < n[i].length; j++) {
//                System.out.println("n(" + i + "," + j + ") = " + n[i][j]);
//            }
//
//        }


        long X = 0;
        // sum X over all k orders

        int[][] korders = new int[C][];
        int[] size = new int[C];

        for (int i = 0; i < C; i++) {
            korders[i] = new int[f[i]];
            size[i] = C - i;
        }

        long startTime = System.currentTimeMillis();

        int calls = 0;
        do {
//            for (int j = 0; j < C; j++) {
//                for (int kj : korders[j]) {
//
//                    System.out.print(kj);
//                }
//                System.out.print("|");
//            }

//            System.out.println();

            for (int i = 0; i < k.length; i++) {
                for (int j = 0; j < k[i].length; j++) {
                    k[i][j] = 0;
                }
            }

            for (int j = 0; j < C; j++) {
                for (int i = 0; i < korders[j].length; i++) {

                    //System.out.println("korders[" + j + ", " + i + "] = " + korders[j][i]);
                    //System.out.println("k(" + korders[j][i] + "," + j + ") += 1");
                    //System.out.println("k.length = " + k.length);

                    k[korders[j][i] + j][j] += 1;
                }
            }

            computeN(k, n, nC, constraints);

            //print(k, n);

            //System.out.println("X = " + X(n, k, f));

            X += X(n, k, f);
            calls += 1;
        }
        while (kton_next(korders, size, 0));

        long stopTime = System.currentTimeMillis();


        System.out.println("Total # constrained ranked histories = " + X + " in " + calls + " calls.");
        System.out.println("Elapsed time " + Math.round((stopTime - startTime) / 10.0) / 100.0 + " seconds");


        //for (int i = 0; i < korders.size(); i++) {
        //
        //    k = korders.get(i);
        //
        //    computeN(k, n, nC, constraints);
        //    X += X(n, k, f);
        //
        //}
        //System.out.println(X + " combinations");


    }
}
