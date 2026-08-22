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
            SericolaRewardDensityWorkspace workspace,
            double[] incSorted,
            double[] out) {

        if (!xIsZero && !xIsOne) {
            cumulants.ensureSecondDifferenceCapacity();
        }

        final double lt = lambda * branchLength;
        final double scaleBase = branchLength * lambda * invAlphaDiff * invAlphaDiff;
        final double[] C = cumulants.values();
        final double[] poissonWeights = workspace.preparePoissonWeights(lt, N);

        for (int n = 0; n <= N; ++n) {
            if (n >= 1) {
                Arrays.fill(incSorted, 0.0);

                if (xIsZero) {
                    addBoundarySecondDifference(h, n, 0, cumulants, C, incSorted);
                } else if (xIsOne) {
                    addBoundarySecondDifference(h, n, n - 1, cumulants, C, incSorted);
                } else {
                    addInteriorSecondDifferences(h, n, xh, cumulants, workspace, incSorted);
                }

                addSortedToOriginal(out, scaleBase * poissonWeights[n] * n, incSorted);
            }
        }
    }

    double contractWrtRewardProportion(
            int h,
            int N,
            double branchLength,
            double lambda,
            double invAlphaDiff,
            double xh,
            boolean xIsZero,
            boolean xIsOne,
            SericolaCumulantMatrices cumulants,
            SericolaRewardDensityWorkspace workspace,
            double[] preOriginal,
            double[] postOriginal) {

        if (preOriginal == null || preOriginal.length != dim) {
            throw new IllegalArgumentException("preOriginal must be length dim=" + dim);
        }
        if (postOriginal == null || postOriginal.length != dim) {
            throw new IllegalArgumentException("postOriginal must be length dim=" + dim);
        }

        if (!xIsZero && !xIsOne) {
            cumulants.ensureSecondDifferenceCapacity();
        }

        final double[] pairWeights = workspace.contractionWeights();
        fillSortedPairWeights(preOriginal, postOriginal, pairWeights);

        final double lt = lambda * branchLength;
        final double scaleBase = branchLength * lambda * invAlphaDiff * invAlphaDiff;
        final double[] C = cumulants.values();
        final double[] poissonWeights = workspace.preparePoissonWeights(lt, N);
        double acc = 0.0;
        double cAcc = 0.0;

        for (int n = 0; n <= N; ++n) {
            if (n >= 1) {
                final double inner;
                if (xIsZero) {
                    inner = contractBoundarySecondDifference(h, n, 0, cumulants, C, pairWeights);
                } else if (xIsOne) {
                    inner = contractBoundarySecondDifference(h, n, n - 1, cumulants, C, pairWeights);
                } else {
                    inner = contractInteriorSecondDifferences(h, n, xh, cumulants, pairWeights, workspace);
                }

                final double y = scaleBase * poissonWeights[n] * n * inner - cAcc;
                final double t = acc + y;
                cAcc = (t - acc) - y;
                acc = t;
            }
        }

        return acc;
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

    private double contractBoundarySecondDifference(
            int h,
            int n,
            int k,
            SericolaCumulantMatrices cumulants,
            double[] C,
            double[] pairWeights) {

        final int c0 = cumulants.offset(h, n + 1, k);
        final int c1 = cumulants.offset(h, n + 1, k + 1);
        final int c2 = cumulants.offset(h, n + 1, k + 2);

        double acc = 0.0;
        int uv = 0;
        for (; uv <= dim2 - 4; uv += 4) {
            acc += pairWeights[uv] * (C[c2 + uv] - 2.0 * C[c1 + uv] + C[c0 + uv]);
            acc += pairWeights[uv + 1] * (C[c2 + uv + 1] - 2.0 * C[c1 + uv + 1] + C[c0 + uv + 1]);
            acc += pairWeights[uv + 2] * (C[c2 + uv + 2] - 2.0 * C[c1 + uv + 2] + C[c0 + uv + 2]);
            acc += pairWeights[uv + 3] * (C[c2 + uv + 3] - 2.0 * C[c1 + uv + 3] + C[c0 + uv + 3]);
        }
        for (; uv < dim2; ++uv) {
            acc += pairWeights[uv] * (C[c2 + uv] - 2.0 * C[c1 + uv] + C[c0 + uv]);
        }
        return acc;
    }

    private void addInteriorSecondDifferences(
            int h,
            int n,
            double xh,
            SericolaCumulantMatrices cumulants,
            SericolaRewardDensityWorkspace workspace,
            double[] incSorted) {

        cumulants.prepareSecondDifferenceRow(h, n);
        final int d2Base = cumulants.secondDifferenceOffset(h, n, 0);
        final double[] d2 = cumulants.secondDifferences();
        final double[] bernsteinWeights = workspace.prepareBernsteinWeights(n - 1, xh);

        for (int k = 0; k <= n - 1; ++k) {
            final int off = d2Base + k * dim2;
            final double wk = bernsteinWeights[k];

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
        }
    }

    private double contractInteriorSecondDifferences(
            int h,
            int n,
            double xh,
            SericolaCumulantMatrices cumulants,
            double[] pairWeights,
            SericolaRewardDensityWorkspace workspace) {

        cumulants.prepareSecondDifferenceRow(h, n);
        final int d2Base = cumulants.secondDifferenceOffset(h, n, 0);
        final double[] d2 = cumulants.secondDifferences();
        final double[] bernsteinWeights = workspace.prepareBernsteinWeights(n - 1, xh);

        double acc = 0.0;
        for (int k = 0; k <= n - 1; ++k) {
            final int off = d2Base + k * dim2;
            final double wk = bernsteinWeights[k];

            double dot = 0.0;
            int uv = 0;
            for (; uv <= dim2 - 4; uv += 4) {
                dot += pairWeights[uv] * d2[off + uv];
                dot += pairWeights[uv + 1] * d2[off + uv + 1];
                dot += pairWeights[uv + 2] * d2[off + uv + 2];
                dot += pairWeights[uv + 3] * d2[off + uv + 3];
            }
            for (; uv < dim2; ++uv) {
                dot += pairWeights[uv] * d2[off + uv];
            }

            acc += wk * dot;
        }

        return acc;
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

    private void fillSortedPairWeights(
            double[] preOriginal,
            double[] postOriginal,
            double[] pairWeights) {

        int uv = 0;
        for (int uS = 0; uS < dim; ++uS) {
            final double pre = preOriginal[outRowBaseBySorted[uS] / dim];
            for (int vS = 0; vS < dim; ++vS) {
                pairWeights[uv++] = pre * postOriginal[outColBySorted[vS]];
            }
        }
    }
}
