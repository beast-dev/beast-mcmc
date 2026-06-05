/*
 * SericolaRewardDensityPdf.java
 *
 * Allocation-free kernel for continuous reward-density accumulation.
 */

package dr.inference.markovjumps;

import java.util.Arrays;

final class SericolaRewardDensityPdf {

    private final int dim;
    private final int dim2;
    private final int[] outRowBaseBySorted;
    private final int[] outColBySorted;

    SericolaRewardDensityPdf(int dim, int[] outRowBaseBySorted, int[] outColBySorted) {
        this.dim = dim;
        this.dim2 = dim * dim;
        this.outRowBaseBySorted = outRowBaseBySorted;
        this.outColBySorted = outColBySorted;
    }

    void accumulateInto(
            double[][] W,
            int T,
            int N,
            double lambda,
            double[] invAlphaDiff,
            SericolaCumulantMatrices cumulants,
            SericolaRewardDensityWorkspace workspace) {

        final double[] C = cumulants.values();
        final int[] H = workspace.intervals();
        final int[] NN = workspace.stepCounts();
        final double[] lt = workspace.lambdaTimes();
        final double[] premult = workspace.premults();
        final double[] ratio = workspace.ratios();
        final double[] w0 = workspace.bernsteinStartWeights();
        final double[] oneMinus = workspace.oneMinus();
        final boolean[] isZero = workspace.isZero();
        final boolean[] isOne = workspace.isOne();
        final double[] incSorted = workspace.increments();

        for (int n = 0; n <= N; ++n) {
            for (int t = 0; t < T; ++t) {
                if (NN[t] < n) {
                    continue;
                }

                final int h = H[t];
                final double time = lt[t] / lambda;
                final double base = (lambda * invAlphaDiff[h]) * premult[t] * time;

                if (isZero[t]) {
                    addDiffBlockToOriginalOrder(
                            W[t],
                            base,
                            cumulants.offset(h, n + 1, 1),
                            cumulants.offset(h, n + 1, 0),
                            C);
                } else if (isOne[t]) {
                    addDiffBlockToOriginalOrder(
                            W[t],
                            base,
                            cumulants.offset(h, n + 1, n + 1),
                            cumulants.offset(h, n + 1, n),
                            C);
                } else {
                    addInteriorToOriginalOrder(
                            h,
                            n,
                            ratio[t],
                            w0[t],
                            premult[t],
                            time,
                            lambda,
                            invAlphaDiff[h],
                            cumulants,
                            C,
                            W[t],
                            incSorted);
                }
            }

            final double inv = 1.0 / (n + 1.0);
            for (int t = 0; t < T; ++t) {
                premult[t] *= lt[t] * inv;
            }

            for (int t = 0; t < T; ++t) {
                if (!isZero[t] && !isOne[t]) {
                    w0[t] *= oneMinus[t];
                }
            }
        }
    }

    private void addInteriorToOriginalOrder(
            int h,
            int n,
            double ratio,
            double w0,
            double premult,
            double time,
            double lambda,
            double invAlphaDiff,
            SericolaCumulantMatrices cumulants,
            double[] C,
            double[] WtOriginal,
            double[] incSorted) {

        final double scale = (lambda * invAlphaDiff) * premult * time;

        Arrays.fill(incSorted, 0.0);

        int aOff = cumulants.offset(h, n + 1, 1);
        int bOff = cumulants.offset(h, n + 1, 0);
        for (int uv = 0; uv < dim2; ++uv) {
            incSorted[uv] += w0 * (C[aOff + uv] - C[bOff + uv]);
        }

        double w = w0;
        for (int k = 0; k < n; ++k) {
            w *= ((double) (n - k) / (double) (k + 1)) * ratio;
            final int kp1 = k + 1;

            aOff = cumulants.offset(h, n + 1, kp1 + 1);
            bOff = cumulants.offset(h, n + 1, kp1);

            for (int uv = 0; uv < dim2; ++uv) {
                incSorted[uv] += w * (C[aOff + uv] - C[bOff + uv]);
            }
        }

        addSortedToOriginal(WtOriginal, scale, incSorted);
    }

    private void addDiffBlockToOriginalOrder(
            double[] WtOriginal,
            double scale,
            int aOff,
            int bOff,
            double[] C) {

        for (int uS = 0; uS < dim; ++uS) {
            final int outRowBase = outRowBaseBySorted[uS];
            final int inRowBase = uS * dim;
            for (int vS = 0; vS < dim; ++vS) {
                final int uvS = inRowBase + vS;
                WtOriginal[outRowBase + outColBySorted[vS]] += scale * (C[aOff + uvS] - C[bOff + uvS]);
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
