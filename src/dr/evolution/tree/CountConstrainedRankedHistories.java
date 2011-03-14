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

                if (constraints[j].intersects(constraints[i])) {
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
        return Rnk;
    }

    private static long X(int[][] n, int[][] k, int[] f) {

        int C = f.length;

        long X = 1;

        // for all levels
        for (int i = 0; i < C; i++) {
            // for all forests in level i
            for (int j = 0; j < i + 1; j++) {
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

            System.out.println("parent = " + parent);

            n[parent][parent] += 1;
        }
    }

    public static void main(String[] args) {

        String[] constraintStrings = {"1111111", "1110000", "0000011"};

        // number of constraints, including the "1,1,...,1,1" constraint
        int C = constraintStrings.length;

        BitSet[] constraints = new BitSet[C];

        for (int i = 0; i < C; i++) {
            constraints[i] = createClade(constraintStrings[i]);
        }

        int[] f = new int[C];

        computeF(constraints, f);

        for (int i = 0; i < f.length; i++) {
            System.out.println((i + 1) + " : " + f[i]);
        }

        System.out.println("Nc = " + Nc(f));

        int[] nC = nC(constraints);
        for (int i = 0; i < nC.length; i++) {
            System.out.println("n_{" + (C - 1) + "," + i + "} = " + nC[i]);
        }


        System.out.println("T(1,2,1) = " + totalOrder(new int[]{1, 2, 1}));

        int[][] k = new int[][]{
                {1}, // 1 node in first level forest
                {0, 0}, // 0 nodes in second level forest of first calibration, 1 nodes in second level forest of second calibration
                {1, 1, 0}
        };

        int[][] n = new int[k.length][];
        for (int i = 0; i < n.length; i++) {
            n[i] = new int[k[i].length];
        }

        computeN(k, n, nC, constraints);

        for (int i = 0; i < n.length; i++) {
            for (int j = 0; j < n[i].length; j++) {
                System.out.println("n(" + i + "," + j + ") = " + n[i][j]);
            }

        }

        List<int[][]> korders = new ArrayList<int[][]>();

        long X = 0;
        // sum X over all k orders
        for (int i = 0; i < korders.size(); i++) {

            k = korders.get(i);

            computeN(k, n, nC, constraints);
            X += X(n, k, f);

        }
        System.out.println(X + " combinations");


    }
}
