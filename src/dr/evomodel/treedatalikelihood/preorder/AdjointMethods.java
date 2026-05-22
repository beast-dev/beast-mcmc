package dr.evomodel.treedatalikelihood.preorder;

import dr.evomodel.substmodel.EigenDecomposition;
import dr.evomodel.substmodel.EigenSystem;
import dr.evomodel.substmodel.LapackEigenSystem;
import dr.math.matrixAlgebra.WrappedVector;

import java.util.Arrays;

public class AdjointMethods {

    public static void main(String[] args) {

        double[][] matrix = {
                {-2, 0, 2},
                {1, -1, 0},
                {0, 1, -1}
        };

        EigenSystem es = new LapackEigenSystem(3);
        EigenDecomposition ed = es.decomposeMatrix(matrix);

        double[] branchLengths = { 0.5 };
        int stateCount = matrix.length;

        double[] transformed = new double[stateCount * stateCount];
        Arrays.fill(transformed, 1.0);

        double[] outRateGradient1 = new double[stateCount * stateCount];
        accumulateEigenBasisGradient1(ed.getEigenValues(), null,
                branchLengths, branchLengths.length, 0,
                transformed, outRateGradient1, stateCount, null);

        System.err.println(new WrappedVector.Raw(outRateGradient1));

        double[] outRateGradient2 = new double[stateCount * stateCount];
        accumulateEigenBasisGradient2(ed.getEigenValues(), null,
                branchLengths, branchLengths.length, 0,
                transformed, outRateGradient2, stateCount, null);

        System.err.println(new WrappedVector.Raw(outRateGradient2));
    }

    static void OneOne(
            double[] eigenValues,
            double t,
            int ls, int rs,
            double[] transformed,
            double[] eigenBasisGrad,
            int S) {

        final double la    = eigenValues[ls];
        final double lbb   = eigenValues[rs];
        final double tla   = t * la;
        final double tlb   = t * lbb;
        final double coeff = (Math.abs(tla - tlb) < 1e-12)
                ? t * Math.exp(tla)
                : (Math.exp(tla) - Math.exp(tlb)) / (la - lbb);

        eigenBasisGrad[ls * S + rs] += transformed[ls * S + rs] * coeff;
    }

    static void OneTwo(
            double[] eigenValues,
            double t,
            int ls, int rs,
            double[] transformed,
            double[] eigenBasisGrad,
            int S) {
        final double la    = eigenValues[ls];
        final double rr    = eigenValues[rs];
        final double ri    = eigenValues[S + rs];
        final double sr    = rr - la;
        final double den   = sr * sr + ri * ri;
        final double scale = Math.exp(t * la);

        double ic0, ic1;
        if (den < 1e-12) {
            ic0 = t;
            ic1 = 0;
        } else {
            final double ex = Math.exp(t * sr);
            final double cs = Math.cos(t * ri);
            final double sn = Math.sin(t * ri);
            ic0 = (ex * (sr * cs + ri * sn) - sr) / den;
            ic1 = (ex * (sr * sn - ri * cs) + ri) / den;
        }

        final double c0  = scale * ic0;
        final double c1  = scale * ic1;
        final double in0 = transformed[ls * S + rs];
        final double in1 = transformed[ls * S + rs + 1];

        eigenBasisGrad[ls * S + rs]     +=  c0 * in0 + c1 * in1;
        eigenBasisGrad[ls * S + rs + 1] += -c1 * in0 + c0 * in1;
    }

    static void TwoOne(
            double[] eigenValues,
            double t,
            int ls, int rs,
            double[] transformed,
            double[] eigenBasisGrad,
            int S) {

        final double lr  = eigenValues[ls];
        final double li  = eigenValues[S + ls];
        final double rb2 = eigenValues[rs];
        final double sr  = rb2 - lr;
        final double den = sr * sr + li * li;

        double ic0, ic1;
        if (den < 1e-12) {
            ic0 = t;
            ic1 = 0;
        } else {
            final double ex = Math.exp(t * sr);
            final double cs = Math.cos(t * li);
            final double sn = Math.sin(t * li);
            ic0 = (ex * (sr * cs + li * sn) - sr) / den;
            ic1 = (ex * (sr * sn - li * cs) + li) / den;
        }

        final double eR = Math.exp(t * lr);
        final double cI = Math.cos(t * li);
        final double sI = Math.sin(t * li);

        final double l00 = eR *  cI;
        final double l01 = eR * -sI;
        final double l10 = eR *  sI;
        final double l11 = eR *  cI;

        final double p0 = l00 * ic0 - l01 * ic1;
        final double p1 = l00 * ic1 + l01 * ic0;
        final double p2 = l10 * ic0 - l11 * ic1;
        final double p3 = l10 * ic1 + l11 * ic0;

        final double in0 = transformed[ls * S + rs];
        final double in1 = transformed[(ls + 1) * S + rs];

        eigenBasisGrad[ls * S + rs]       += p0 * in0 + p1 * in1;
        eigenBasisGrad[(ls + 1) * S + rs] += p2 * in0 + p3 * in1;
    }

    static void TwoTwo(
            double[] eigenValues,
            double t,
            int ls, int rs,
            double[] transformed,
            double[] eigenBasisGrad,
            int S) {

        final double lr = eigenValues[ls];
        final double li = eigenValues[S + ls];
        final double rr = eigenValues[rs];
        final double ri = eigenValues[S + rs];

        final double sr1 = rr - lr;
        final double si1 = li + ri;
        final double sr2 = rr - lr;
        final double si2 = ri - li;

        final double e1    = Math.exp(t * lr);
        final double exp1r = e1 * Math.cos(-t * li);
        final double exp1i = e1 * Math.sin(-t * li);
        final double exp2r = e1 * Math.cos( t * li);
        final double exp2i = e1 * Math.sin( t * li);

        double int1r, int1i;
        final double d1 = sr1 * sr1 + si1 * si1;
        if (d1 < 1e-12) {
            int1r = t;
            int1i = 0;
        } else {
            final double ex1 = Math.exp(t * sr1);
            final double cs1 = Math.cos(t * si1);
            final double sn1 = Math.sin(t * si1);
            int1r = (sr1 * (ex1 * cs1 - 1) + si1 * ex1 * sn1) / d1;
            int1i = (sr1 *  ex1 * sn1 - si1 * (ex1 * cs1 - 1)) / d1;
        }

        double int2r, int2i;
        final double d2 = sr2 * sr2 + si2 * si2;
        if (d2 < 1e-12) {
            int2r = t;
            int2i = 0;
        } else {
            final double ex2 = Math.exp(t * sr2);
            final double cs2 = Math.cos(t * si2);
            final double sn2 = Math.sin(t * si2);
            int2r = (sr2 * (ex2 * cs2 - 1) + si2 * ex2 * sn2) / d2;
            int2i = (sr2 *  ex2 * sn2 - si2 * (ex2 * cs2 - 1)) / d2;
        }

        final double pr  = exp1r * int1r - exp1i * int1i;
        final double pi_ = exp1r * int1i + exp1i * int1r;
        final double mr_ = exp2r * int2r - exp2i * int2i;
        final double mi_ = exp2r * int2i + exp2i * int2r;

        final double[][] basis = {
                { .5,  0,  .5,  0},
                {  0, .5,   0, .5},
                {  0,-.5,   0, .5},
                { .5,  0, -.5,  0},
        };

        double[] c = new double[16];
        for (int col = 0; col < 4; col++) {
            final double u = basis[col][0];
            final double v = basis[col][1];
            final double p = basis[col][2];
            final double q = basis[col][3];
            c[     col] =  mr_ * u + mi_ * v + pr  * p + pi_ * q;
            c[ 4 + col] = -mi_ * u + mr_ * v - pi_ * p + pr  * q;
            c[ 8 + col] =  mi_ * u - mr_ * v - pi_ * p + pr  * q;
            c[12 + col] =  mr_ * u + mi_ * v - pr  * p - pi_ * q;
        }

        final double in00 = transformed[ls * S + rs];
        final double in01 = transformed[ls * S + rs + 1];
        final double in10 = transformed[(ls + 1) * S + rs];
        final double in11 = transformed[(ls + 1) * S + rs + 1];

        eigenBasisGrad[ls * S + rs]           += c[ 0]*in00 + c[ 1]*in01 + c[ 2]*in10 + c[ 3]*in11;
        eigenBasisGrad[ls * S + rs + 1]       += c[ 4]*in00 + c[ 5]*in01 + c[ 6]*in10 + c[ 7]*in11;
        eigenBasisGrad[(ls + 1) * S + rs]     += c[ 8]*in00 + c[ 9]*in01 + c[10]*in10 + c[11]*in11;
        eigenBasisGrad[(ls + 1) * S + rs + 1] += c[12]*in00 + c[13]*in01 + c[14]*in10 + c[15]*in11;
    }

    static int accumulateEigenBasisGradient1(
            double[] eigenValues,
            int[]    matrixIndices,
            double[] branchLengths,
            int      matrixCount,
            int      hasComplexEigenvalues,
            double[] transformed,
            double[] outRateGradient,
            int      kStateCount,
            double[][] gTransitionMatrices) {

        final int S  = kStateCount;
        final int S2 = S * S;
        final int M  = matrixCount;

        int[] blockStarts  = new int[S];
        int[] blockDimsArr = new int[S];

        int numBlocks = 0;
        for (int i = 0; i < S; i++) {
            blockStarts[numBlocks] = i;
            if (Math.abs(eigenValues[S + i]) > 1e-12) {
                blockDimsArr[numBlocks] = 2;
                i++;
            } else {
                blockDimsArr[numBlocks] = 1;
            }
            numBlocks++;
        }

        double[] eigenBasisGrad = outRateGradient;
        Arrays.fill(eigenBasisGrad, 0, S2, 0.0);

        for (int m = 0; m < M; m++) {

            final double t = branchLengths[m];

            for (int lb = 0; lb < numBlocks; lb++) {
                final int ls = blockStarts[lb], ld = blockDimsArr[lb];

                for (int rb = 0; rb < numBlocks; rb++) {
                    final int rs = blockStarts[rb], rd = blockDimsArr[rb];

                    if (ld == 1 && rd == 1) {

                        OneOne(eigenValues, t, ls, rs, transformed, eigenBasisGrad, S);

                    } else if (ld == 1 && rd == 2) {

                        OneTwo(eigenValues, t, ls, rs, transformed, eigenBasisGrad, S);

                    } else if (ld == 2 && rd == 1) {

                        TwoOne(eigenValues, t, ls, rs, transformed, eigenBasisGrad, S);

                    } else {

                        TwoTwo(eigenValues, t, ls, rs, transformed, eigenBasisGrad, S);

                    }
                }
            }
        }
        return 0;
    }

    static int accumulateEigenBasisGradient2(
            double[] eigenValues,
            int[]    matrixIndices,
            double[] branchLengths,
            int      matrixCount,
            int      hasComplexEigenvalues,
            double[] transformed,
            double[] outRateGradient,
            int      kStateCount,
            double[][] gTransitionMatrices) {

        final int S  = kStateCount;
        final int S2 = S * S;
        final int M  = matrixCount;

        double[] eigenBasisGrad = outRateGradient;
        Arrays.fill(eigenBasisGrad, 0, S2, 0.0);

        for (int m = 0; m < M; m++) {
            final double t = branchLengths[m];

            for (int ls = 0; ls < S; ++ls) {
                if (eigenValues[S + ls] == 0.0) {
                    for (int rs = 0; rs < S; ++rs) {
                        if (eigenValues[S + rs] == 0.0) {

                            OneOne(eigenValues, t, ls, rs, transformed, eigenBasisGrad, S);

                        } else {

                            OneTwo(eigenValues, t, ls, rs, transformed, eigenBasisGrad, S);

                            ++rs;
                        }
                    }
                } else {
                    for (int rs = 0; rs < S; ++rs) {
                        if (eigenValues[S + rs] == 0.0) {

                            TwoOne(eigenValues, t, ls, rs, transformed, eigenBasisGrad, S);

                        } else {

                            TwoTwo(eigenValues, t, ls, rs, transformed, eigenBasisGrad, S);

                            ++rs;
                        }
                    }

                    ++ls;
                }
            }

        }
        return 0;
    }
}
