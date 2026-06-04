/*
 * SericolaRewardDensityDerivative.java
 *
 * Allocation-free kernel for the reward-density derivative with respect to rewardProportion.
 */

package dr.inference.markovjumps;

import java.util.Arrays;

final class SericolaRewardDensityDerivative {

    private final int dim;
    private final int dim2;
    private final int[] outRowBaseBySorted;
    private final int[] outColBySorted;

    SericolaRewardDensityDerivative(int dim, int[] outRowBaseBySorted, int[] outColBySorted) {
        this.dim = dim;
        this.dim2 = dim * dim;
        this.outRowBaseBySorted = outRowBaseBySorted;
        this.outColBySorted = outColBySorted;
    }

    void computeWrtRewardProportionInto(
            int h,
            int N,
            double branchLength,
            double lambda,
            double invAlphaDiff,
            double xh,
            boolean xIsZero,
            boolean xIsOne,
            SericolaCumulantMatrices cumulants,
            double[] incSorted,
            double[] out) {

        if (!xIsZero && !xIsOne) {
            cumulants.ensureSecondDifferenceCapacity();
        }

        final double lt = lambda * branchLength;
        final double oneMinus = 1.0 - xh;
        final double ratio = (!xIsZero && !xIsOne) ? (xh / oneMinus) : 0.0;
        final double scaleBase = branchLength * lambda * invAlphaDiff * invAlphaDiff;
        final double[] C = cumulants.values();

        double premult = Math.exp(-lt);
        double w0m = 1.0;

        for (int n = 0; n <= N; ++n) {
            if (n >= 1) {
                Arrays.fill(incSorted, 0.0);

                if (xIsZero) {
                    addBoundarySecondDifference(h, n, 0, cumulants, C, incSorted);
                } else if (xIsOne) {
                    addBoundarySecondDifference(h, n, n - 1, cumulants, C, incSorted);
                } else {
                    addInteriorSecondDifferences(h, n, ratio, w0m, cumulants, incSorted);
                }

                addSortedToOriginal(out, scaleBase * premult * n, incSorted);
            }

            premult *= lt / (n + 1.0);
            if (n >= 1) {
                w0m *= oneMinus;
            }
        }
    }

    private void addBoundarySecondDifference(
            int h,
            int n,
            int k,
            SericolaCumulantMatrices cumulants,
            double[] C,
            double[] incSorted) {

        final int c0 = cumulants.offset(h, n + 1, k);
        final int c1 = cumulants.offset(h, n + 1, k + 1);
        final int c2 = cumulants.offset(h, n + 1, k + 2);

        for (int uv = 0; uv < dim2; ++uv) {
            incSorted[uv] += C[c2 + uv] - 2.0 * C[c1 + uv] + C[c0 + uv];
        }
    }

    private void addInteriorSecondDifferences(
            int h,
            int n,
            double ratio,
            double w0m,
            SericolaCumulantMatrices cumulants,
            double[] incSorted) {

        cumulants.prepareSecondDifferenceRow(h, n);
        final int d2Base = cumulants.secondDifferenceOffset(h, n, 0);
        final double[] d2 = cumulants.secondDifferences();

        double w = w0m;
        for (int k = 0; k <= n - 1; ++k) {
            final int off = d2Base + k * dim2;
            final double wk = w;

            int uv = 0;
            for (; uv <= dim2 - 4; uv += 4) {
                incSorted[uv] += wk * d2[off + uv];
                incSorted[uv + 1] += wk * d2[off + uv + 1];
                incSorted[uv + 2] += wk * d2[off + uv + 2];
                incSorted[uv + 3] += wk * d2[off + uv + 3];
            }
            for (; uv < dim2; ++uv) {
                incSorted[uv] += wk * d2[off + uv];
            }

            if (k < n - 1) {
                w *= ((double) (n - 1 - k) / (double) (k + 1)) * ratio;
            }
        }
    }

    private void addSortedToOriginal(double[] out, double scale, double[] incSorted) {
        for (int uS = 0; uS < dim; ++uS) {
            final int outRowBase = outRowBaseBySorted[uS];
            final int inRowBase = uS * dim;
            for (int vS = 0; vS < dim; ++vS) {
                out[outRowBase + outColBySorted[vS]] += scale * incSorted[inRowBase + vS];
            }
        }
    }
}
