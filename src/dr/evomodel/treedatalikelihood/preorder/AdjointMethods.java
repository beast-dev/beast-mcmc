package dr.evomodel.treedatalikelihood.preorder;

import dr.evomodel.substmodel.EigenDecomposition;
import dr.evomodel.substmodel.EigenSystem;
import dr.evomodel.substmodel.LapackEigenSystem;
import dr.math.MathUtils;
import dr.math.matrixAlgebra.WrappedVector;

import java.util.Arrays;

public class AdjointMethods {

    public static void main(String[] args) {

        MathUtils.setSeed(666);

        double[][] matrix = {
                {-2, 0, 2},
                {1, -1, 0},
                {0, 1, -1}
        };

        EigenSystem es = new LapackEigenSystem(3);
        EigenDecomposition ed = es.decomposeMatrix(matrix);

        double[] branchLengths = { 0.5 };
        int stateCount = matrix.length;
        double scale = 2.0;

        double[] lhs = makeRandomVector(stateCount);
        double[] rhs = makeRandomVector(stateCount);
        double[] transformed = outerProduct(lhs, rhs, scale, stateCount);

//        double[] transformed = new double[stateCount * stateCount];
//        Arrays.fill(transformed, 1.0);
//
//        double[] lhs = new double[stateCount];
//        Arrays.fill(lhs, 1.0);
//        double[] rhs = new double[stateCount];
//        Arrays.fill(rhs, 1.0);
//        double scale = 1.0;

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

        double[] outRateGradient3 = new double[stateCount * stateCount];
        accumulateEigenBasisGradient3(ed.getEigenValues(), null,
                branchLengths, branchLengths.length, 0,
                lhs, rhs, scale,
                outRateGradient3, stateCount, null);

        System.err.println(new WrappedVector.Raw(outRateGradient3));
    }

    static double[] makeRandomVector(int dim) {
        double[] vec = new double[dim];
        for (int i = 0; i < dim; ++i) {
            vec[i] = MathUtils.nextGaussian();
        }
        return vec;
    }

    static double[] outerProduct(double[] lhs, double[] rhs, double scale, int stride) {
        final int dim = lhs.length;

        double[] mat = new double[dim * stride];
        for (int i = 0; i < dim; ++i) {
            for (int j = 0; j < dim; ++j) {
                mat[i * stride + j] = scale * lhs[i] * rhs[j];
            }
        }
        return mat;
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

    static void OneOneNew(
            double[] eigenValues,
            double t,
            int ls, int rs,
//            double[] transformed,
            double lv,
            double rv,
            double scale,
            double[] eigenBasisGrad,
            int S) {

        final double la    = eigenValues[ls];
        final double lbb   = eigenValues[rs];
        final double tla   = t * la;
        final double tlb   = t * lbb;
        final double coeff = (Math.abs(tla - tlb) < 1e-12)
                ? t * Math.exp(tla)
                : (Math.exp(tla) - Math.exp(tlb)) / (la - lbb);

        eigenBasisGrad[ls * S + rs] += lv * rv * scale * coeff;
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

    // OneTwoNew: no duplicate computation exists in OneTwo (the left eigenvalue is real, so its
    // scale factor exp(t*la) involves no trig). Structurally identical to OneTwo.
    static void OneTwoNew(
            double[] eigenValues,
            double t,
            int ls, int rs,
            double lv1, double rv1, double rv2,
            double sc,
//            double[] transformed,
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
        final double in0 = lv1 * rv1 * sc; //transformed[ls * S + rs];
        final double in1 = lv1 * rv2 * sc; //transformed[ls * S + rs + 1];

        eigenBasisGrad[ls * S + rs]     +=  c0 * in0 + c1 * in1;
        eigenBasisGrad[ls * S + rs + 1] += -c1 * in0 + c0 * in1;
    }

    // TwoOneNew: hoists cos/sin of (t*li) before the if/else, eliminating the duplicate
    // computation that exists in TwoOne (lines 118-119 vs 125-126 in the original).
    static void TwoOneNew(
            double[] eigenValues,
            double t,
            int ls, int rs,
//            double[] transformed,
            double lv1, double lv2, double rv1, double sc,
            double[] eigenBasisGrad,
            int S) {

        final double lr  = eigenValues[ls];
        final double li  = eigenValues[S + ls];
        final double rb2 = eigenValues[rs];
        final double sr  = rb2 - lr;
        final double den = sr * sr + li * li;

        final double eR = Math.exp(t * lr);
        final double cI = Math.cos(t * li);
        final double sI = Math.sin(t * li);

        double ic0, ic1;
        if (den < 1e-12) {
            ic0 = t;
            ic1 = 0;
        } else {
            final double ex = Math.exp(t * sr);
            ic0 = (ex * (sr * cI + li * sI) - sr) / den;
            ic1 = (ex * (sr * sI - li * cI) + li) / den;
        }

        final double l00 = eR *  cI;
        final double l01 = eR * -sI;
        final double l10 = eR *  sI;
        final double l11 = eR *  cI;

        final double p0 = l00 * ic0 - l01 * ic1;
        final double p1 = l00 * ic1 + l01 * ic0;
        final double p2 = l10 * ic0 - l11 * ic1;
        final double p3 = l10 * ic1 + l11 * ic0;

        final double in0 = lv1 * rv1 * sc; // transformed[ls * S + rs];
        final double in1 = lv2 * rv1 * sc; // transformed[(ls + 1) * S + rs];

        eigenBasisGrad[ls * S + rs]       += p0 * in0 + p1 * in1;
        eigenBasisGrad[(ls + 1) * S + rs] += p2 * in0 + p3 * in1;
    }

    // TwoTwoNew applies four optimizations over TwoTwo:
    //   1. sr1 == sr2 always (both = rr - lr); merged into a single sr.
    //   2. cos(-t*li) == cos(t*li) and sin(-t*li) == -sin(t*li), so exp1r == exp2r and
    //      exp1i == -exp2i; compute er = e1*cos(t*li) and ei = e1*sin(t*li) once.
    //   3. ex = exp(t*sr) is shared between the two integral blocks; computed once.
    //   4. The basis[][] and c[] heap allocations are replaced by four scalars A,B,C,D
    //      derived by unrolling the constant basis matrix.
    static void TwoTwoNew(
            double[] eigenValues,
            double t,
            int ls, int rs,
//            double[] transformed,
            double lv1, double lv2,
            double rv1, double rv2,
            double sc,
            double[] eigenBasisGrad,
            int S) {

        final double lr = eigenValues[ls];
        final double li = eigenValues[S + ls];
        final double rr = eigenValues[rs];
        final double ri = eigenValues[S + rs];

        final double sr  = rr - lr;
        final double si1 = li + ri;
        final double si2 = ri - li;

        final double e1 = Math.exp(t * lr);
        final double er = e1 * Math.cos(t * li);   // exp1r = exp2r  (cos is even)
        final double ei = e1 * Math.sin(t * li);   // exp2i;  exp1i = -ei

        final double sr2 = sr * sr;
        final double d1  = sr2 + si1 * si1;
        final double d2  = sr2 + si2 * si2;

        final double ex = (d1 >= 1e-12 || d2 >= 1e-12) ? Math.exp(t * sr) : 0.0;

        double int1r, int1i;
        if (d1 < 1e-12) {
            int1r = t;
            int1i = 0;
        } else {
            final double cs1       = Math.cos(t * si1);
            final double sn1       = Math.sin(t * si1);
            final double ex_cs1_m1 = ex * cs1 - 1;
            final double ex_sn1    = ex * sn1;
            int1r = (sr * ex_cs1_m1 + si1 * ex_sn1) / d1;
            int1i = (sr * ex_sn1    - si1 * ex_cs1_m1) / d1;
        }

        double int2r, int2i;
        if (d2 < 1e-12) {
            int2r = t;
            int2i = 0;
        } else {
            final double cs2       = Math.cos(t * si2);
            final double sn2       = Math.sin(t * si2);
            final double ex_cs2_m1 = ex * cs2 - 1;
            final double ex_sn2    = ex * sn2;
            int2r = (sr * ex_cs2_m1 + si2 * ex_sn2) / d2;
            int2i = (sr * ex_sn2    - si2 * ex_cs2_m1) / d2;
        }

        // exp1i = -ei, so pr/pi_ simplify vs the original
        final double pr  = er * int1r + ei * int1i;
        final double pi_ = er * int1i - ei * int1r;
        final double mr_ = er * int2r - ei * int2i;
        final double mi_ = er * int2i + ei * int2r;

        // Unrolled product of the constant basis matrix with [mr_, mi_, pr, pi_]:
        //   A = 0.5*(mr_+pr), B = 0.5*(mi_+pi_), C = 0.5*(pi_-mi_), D = 0.5*(mr_-pr)
        // giving the 4x4 output matrix [ A  B  C  D / -B  A -D  C / -C -D  A  B / D -C -B  A ]
        final double A = 0.5 * (mr_ + pr);
        final double B = 0.5 * (mi_ + pi_);
        final double C = 0.5 * (pi_ - mi_);
        final double D = 0.5 * (mr_ - pr);

        final double in00 = lv1 * rv1 * sc; // transformed[ls * S + rs];
        final double in01 = lv1 * rv2 * sc; // transformed[ls * S + rs + 1];
        final double in10 = lv2 * rv1 * sc; // transformed[(ls + 1) * S + rs];
        final double in11 = lv2 * rv2 * sc; // transformed[(ls + 1) * S + rs + 1];

        eigenBasisGrad[ls * S + rs]           +=  A*in00 + B*in01 + C*in10 + D*in11;
        eigenBasisGrad[ls * S + rs + 1]       += -B*in00 + A*in01 - D*in10 + C*in11;
        eigenBasisGrad[(ls + 1) * S + rs]     += -C*in00 - D*in01 + A*in10 + B*in11;
        eigenBasisGrad[(ls + 1) * S + rs + 1] +=  D*in00 - C*in01 - B*in10 + A*in11;
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

    static int accumulateEigenBasisGradient3(
            double[] eigenValues,
            int[]    matrixIndices,
            double[] branchLengths,
            int      matrixCount,
            int      hasComplexEigenvalues,
//            double[] transformed,
            double[] lhs,
            double[] rhs,
            double scale,
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
                final double lv1 = lhs[ls];

                if (eigenValues[S + ls] == 0.0) {
                    for (int rs = 0; rs < S; ++rs) {
                        final double rv1 = rhs[rs];

                        if (eigenValues[S + rs] == 0.0) {

                            OneOneNew(eigenValues, t, ls, rs, lv1, rv1, scale, eigenBasisGrad, S);

                        } else {

                            final double rv2 = rhs[rs + 1];
                            OneTwoNew(eigenValues, t, ls, rs, lv1, rv1, rv2, scale, eigenBasisGrad, S);

                            ++rs;
                        }
                    }
                } else {
                    final double lv2 = lhs[ls + 1];

                    for (int rs = 0; rs < S; ++rs) {
                        final double rv1 = rhs[rs];

                        if (eigenValues[S + rs] == 0.0) {

                            TwoOneNew(eigenValues, t, ls, rs, lv1, lv2, rv1, scale, eigenBasisGrad, S);

                        } else {

                            final double rv2 = rhs[rs + 1];
                            TwoTwoNew(eigenValues, t, ls, rs, lv1, lv2, rv1, rv2, scale, eigenBasisGrad, S);

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
