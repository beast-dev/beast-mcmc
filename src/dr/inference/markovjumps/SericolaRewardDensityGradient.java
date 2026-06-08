/*
 * SericolaRewardDensityGradient.java
 *
 * Reverse-mode kernel for reward-density gradients with respect to the
 * uniformization transition matrix and reward-rate proportions.
 */

package dr.inference.markovjumps;

import java.util.Arrays;

final class SericolaRewardDensityGradient {

    private final int dim;
    private final int dim2;
    private final int phi;
    private final int[] outRowBaseBySorted;
    private final int[] outColBySorted;

    private double[] cumulantAdjoints;
    private double[] transitionMatrixAdjoint;
    private double[] transitionPowers;
    private double[] transitionPowerAdjoints;
    private double[] sortedDensityAdjoint;
    private double[] rewardRateAdjoint;
    private int allocatedN = -1;

    SericolaRewardDensityGradient(int dim, int[] outRowBaseBySorted, int[] outColBySorted) {
        this.dim = dim;
        this.dim2 = dim * dim;
        this.phi = dim - 1;
        this.outRowBaseBySorted = outRowBaseBySorted;
        this.outColBySorted = outColBySorted;
    }

    double computeWrtUniformizationMatrixInto(
            double[] densityAdjointOriginal,
            int h,
            int N,
            double time,
            double lambda,
            double invAlphaDiff,
            boolean xIsZero,
            boolean xIsOne,
            double xh,
            boolean conditionalOnZ0,
            double[] transitionMatrix,
            double[] sortedAlpha,
            SericolaCumulantMatrices cumulants) {

        if (densityAdjointOriginal == null || densityAdjointOriginal.length != dim2) {
            throw new IllegalArgumentException("densityAdjointOriginal must be length dim*dim=" + dim2);
        }
        if (N < 0) {
            throw new IllegalArgumentException("N must be >= 0");
        }

        final int maxExternal = N + 1;
        ensureCapacity(cumulants.allocatedN());
        clearAdjoints(maxExternal);

        final double directLambdaAdjoint = seedCumulantAdjointsFromDensity(
                densityAdjointOriginal,
                h,
                N,
                time,
                lambda,
                invAlphaDiff,
                xIsZero,
                xIsOne,
                xh,
                conditionalOnZ0,
                false,
                sortedAlpha,
                cumulants);

        computeTransitionPowers(maxExternal, transitionMatrix);
        reverseCumulants(maxExternal, transitionMatrix, sortedAlpha, cumulants, true, false);
        reverseTransitionPowers(maxExternal, transitionMatrix);

        return directLambdaAdjoint;
    }

    void computeWrtRewardRatesInto(
            double[] densityAdjointOriginal,
            int h,
            int N,
            double time,
            double lambda,
            double invAlphaDiff,
            boolean xIsZero,
            boolean xIsOne,
            double xh,
            boolean conditionalOnZ0,
            double[] transitionMatrix,
            double[] sortedAlpha,
            SericolaCumulantMatrices cumulants,
            double[] outSortedRewardRateAdjoint) {

        if (densityAdjointOriginal == null || densityAdjointOriginal.length != dim2) {
            throw new IllegalArgumentException("densityAdjointOriginal must be length dim*dim=" + dim2);
        }
        if (outSortedRewardRateAdjoint == null || outSortedRewardRateAdjoint.length != dim) {
            throw new IllegalArgumentException("outSortedRewardRateAdjoint must be length dim=" + dim);
        }
        if (N < 0) {
            throw new IllegalArgumentException("N must be >= 0");
        }
        if (xIsZero || xIsOne) {
            throw new UnsupportedOperationException(
                    "Reward-rate gradients at reward-rate breakpoints are not defined by the interior " +
                            "reverse-mode convention. xh=" + xh + ", h=" + h);
        }

        final int maxExternal = N + 1;
        ensureCapacity(cumulants.allocatedN());
        clearAdjoints(maxExternal);

        seedCumulantAdjointsFromDensity(
                densityAdjointOriginal,
                h,
                N,
                time,
                lambda,
                invAlphaDiff,
                xIsZero,
                xIsOne,
                xh,
                conditionalOnZ0,
                true,
                sortedAlpha,
                cumulants);

        reverseCumulants(maxExternal, transitionMatrix, sortedAlpha, cumulants, false, true);
        System.arraycopy(rewardRateAdjoint, 0, outSortedRewardRateAdjoint, 0, dim);
    }

    double[] transitionMatrixAdjoint() {
        return transitionMatrixAdjoint;
    }

    private double seedCumulantAdjointsFromDensity(
            double[] densityAdjointOriginal,
            int h,
            int N,
            double time,
            double lambda,
            double invAlphaDiff,
            boolean xIsZero,
            boolean xIsOne,
            double xh,
            boolean conditionalOnZ0,
            boolean collectRewardRateAdjoints,
            double[] sortedAlpha,
            SericolaCumulantMatrices cumulants) {

        final double lambdaTime = lambda * time;
        final double denominator;
        final double rawAdjointScale;
        if (conditionalOnZ0) {
            denominator = -Math.expm1(-lambdaTime);
            if (denominator <= 0.0) {
                throw new IllegalStateException("Invalid denominator for conditional probability: " + denominator);
            }
            rawAdjointScale = 1.0 / denominator;
        } else {
            denominator = 1.0;
            rawAdjointScale = 1.0;
        }

        mapDensityAdjointToSorted(densityAdjointOriginal, rawAdjointScale);

        final double[] C = cumulants.values();
        double premult = Math.exp(-lambdaTime);
        double previousPremult = 0.0;
        double directLambdaAdjoint = 0.0;
        double conditionalObjective = 0.0;

        final double oneMinus = 1.0 - xh;
        final double ratio = (!xIsZero && !xIsOne) ? (xh / oneMinus) : 0.0;
        double w0 = 1.0;

        for (int n = 0; n <= N; ++n) {
            final double scale = lambda * invAlphaDiff * premult * time;
            final double scaleLambdaDerivative =
                    invAlphaDiff * time * (premult + lambdaTime * (previousPremult - premult));
            final double innerDot;

            if (xIsZero) {
                innerDot = addDiffCumulantAdjoints(
                        scale,
                        cumulants.offset(h, n + 1, 1),
                        cumulants.offset(h, n + 1, 0),
                        C);
            } else if (xIsOne) {
                innerDot = addDiffCumulantAdjoints(
                        scale,
                        cumulants.offset(h, n + 1, n + 1),
                        cumulants.offset(h, n + 1, n),
                        C);
            } else {
                innerDot = addInteriorCumulantAdjoints(h, n, ratio, w0, scale, cumulants, C);
            }

            directLambdaAdjoint += scaleLambdaDerivative * innerDot;
            if (conditionalOnZ0) {
                conditionalObjective += scale * innerDot;
            }
            if (collectRewardRateAdjoints) {
                final double intervalScaleAdjoint = scale * invAlphaDiff * innerDot;
                rewardRateAdjoint[h - 1] += intervalScaleAdjoint;
                rewardRateAdjoint[h] -= intervalScaleAdjoint;
            }

            previousPremult = premult;
            premult *= lambdaTime / (n + 1.0);
            if (!xIsZero && !xIsOne) {
                w0 *= oneMinus;
            }
        }

        if (conditionalOnZ0) {
            directLambdaAdjoint -= conditionalObjective * time * Math.exp(-lambdaTime) / denominator;
        }

        if (collectRewardRateAdjoints) {
            final double xhAdjoint = computeXhAdjointFromDensity(
                    h,
                    N,
                    time,
                    lambda,
                    invAlphaDiff,
                    xh,
                    cumulants);
            rewardRateAdjoint[h - 1] += xhAdjoint * (xh - 1.0) * invAlphaDiff;
            rewardRateAdjoint[h] -= xhAdjoint * xh * invAlphaDiff;
        }

        return directLambdaAdjoint;
    }

    private void mapDensityAdjointToSorted(double[] densityAdjointOriginal, double scale) {
        for (int uS = 0; uS < dim; ++uS) {
            final int originalRowBase = outRowBaseBySorted[uS];
            final int sortedRowBase = uS * dim;
            for (int vS = 0; vS < dim; ++vS) {
                sortedDensityAdjoint[sortedRowBase + vS] =
                        scale * densityAdjointOriginal[originalRowBase + outColBySorted[vS]];
            }
        }
    }

    private double addInteriorCumulantAdjoints(
            int h,
            int n,
            double ratio,
            double w0,
            double scale,
            SericolaCumulantMatrices cumulants,
            double[] C) {

        double innerDot = 0.0;
        double w = w0;

        for (int k = 0; k <= n; ++k) {
            innerDot += w * addDiffCumulantAdjoints(
                    scale * w,
                    cumulants.offset(h, n + 1, k + 1),
                    cumulants.offset(h, n + 1, k),
                    C);

            if (k < n) {
                w *= ((double) (n - k) / (double) (k + 1)) * ratio;
            }
        }

        return innerDot;
    }

    private double addDiffCumulantAdjoints(double scale, int plusOffset, int minusOffset, double[] C) {
        double innerDot = 0.0;

        for (int uv = 0; uv < dim2; ++uv) {
            final double adjoint = sortedDensityAdjoint[uv];
            final double scaledAdjoint = scale * adjoint;
            cumulantAdjoints[plusOffset + uv] += scaledAdjoint;
            cumulantAdjoints[minusOffset + uv] -= scaledAdjoint;
            innerDot += adjoint * (C[plusOffset + uv] - C[minusOffset + uv]);
        }

        return innerDot;
    }

    private double computeXhAdjointFromDensity(
            int h,
            int N,
            double time,
            double lambda,
            double invAlphaDiff,
            double xh,
            SericolaCumulantMatrices cumulants) {

        cumulants.ensureSecondDifferenceCapacity();

        final double lambdaTime = lambda * time;
        final double oneMinus = 1.0 - xh;
        final double ratio = xh / oneMinus;
        final double scaleBase = time * lambda * invAlphaDiff;

        double premult = Math.exp(-lambdaTime);
        double w0m = 1.0;
        double xhAdjoint = 0.0;

        for (int n = 0; n <= N; ++n) {
            if (n >= 1) {
                final double innerDot = secondDifferenceInnerDot(h, n, ratio, w0m, cumulants);
                xhAdjoint += scaleBase * premult * n * innerDot;
            }

            premult *= lambdaTime / (n + 1.0);
            if (n >= 1) {
                w0m *= oneMinus;
            }
        }

        return xhAdjoint;
    }

    private double secondDifferenceInnerDot(
            int h,
            int n,
            double ratio,
            double w0m,
            SericolaCumulantMatrices cumulants) {

        cumulants.prepareSecondDifferenceRow(h, n);
        final double[] d2 = cumulants.secondDifferences();
        final int baseOffset = cumulants.secondDifferenceOffset(h, n, 0);

        double innerDot = 0.0;
        double w = w0m;

        for (int k = 0; k <= n - 1; ++k) {
            final int off = baseOffset + k * dim2;
            double entryDot = 0.0;

            for (int uv = 0; uv < dim2; ++uv) {
                entryDot += sortedDensityAdjoint[uv] * d2[off + uv];
            }
            innerDot += w * entryDot;

            if (k < n - 1) {
                w *= ((double) (n - 1 - k) / (double) (k + 1)) * ratio;
            }
        }

        return innerDot;
    }

    private void reverseCumulants(
            int maxExternal,
            double[] transitionMatrix,
            double[] sortedAlpha,
            SericolaCumulantMatrices cumulants,
            boolean collectTransitionAdjoints,
            boolean collectRewardRateAdjoints) {

        for (int n = maxExternal; n >= 1; --n) {
            reverseBackwardSweep(n, transitionMatrix, sortedAlpha, cumulants,
                    collectTransitionAdjoints, collectRewardRateAdjoints);
            reverseTerminalOverwrite(n, cumulants, collectTransitionAdjoints);
            reverseForwardSweep(n, transitionMatrix, sortedAlpha, cumulants,
                    collectTransitionAdjoints, collectRewardRateAdjoints);
        }
    }

    private void reverseBackwardSweep(
            int n,
            double[] transitionMatrix,
            double[] sortedAlpha,
            SericolaCumulantMatrices cumulants,
            boolean collectTransitionAdjoints,
            boolean collectRewardRateAdjoints) {

        for (int h = 1; h <= phi; ++h) {
            for (int k = 0; k <= n - 1; ++k) {
                reverseLowRows(h, n, k, transitionMatrix, sortedAlpha, cumulants,
                        collectTransitionAdjoints, collectRewardRateAdjoints);
            }

            if (h < phi) {
                reverseLowBoundaryCopy(h, n, cumulants);
            }
        }
    }

    private void reverseForwardSweep(
            int n,
            double[] transitionMatrix,
            double[] sortedAlpha,
            SericolaCumulantMatrices cumulants,
            boolean collectTransitionAdjoints,
            boolean collectRewardRateAdjoints) {

        for (int h = phi; h >= 1; --h) {
            if (h < phi) {
                reverseHighBoundaryCopy(h, n, cumulants);
            }

            for (int k = n; k >= 1; --k) {
                reverseHighRows(h, n, k, transitionMatrix, sortedAlpha, cumulants,
                        collectTransitionAdjoints, collectRewardRateAdjoints);
            }
        }
    }

    private void reverseLowRows(
            int h,
            int n,
            int k,
            double[] transitionMatrix,
            double[] sortedAlpha,
            SericolaCumulantMatrices cumulants,
            boolean collectTransitionAdjoints,
            boolean collectRewardRateAdjoints) {

        final double[] C = cumulants.values();
        final int curOff = cumulants.offset(h, n, k);
        final int curOffKp1 = cumulants.offset(h, n, k + 1);
        final int prevOff = cumulants.offset(h, n - 1, k);

        for (int u = 0; u <= h - 1; ++u) {
            final int uRow = u * dim;
            final double aScalar = (sortedAlpha[h - 1] - sortedAlpha[u]) /
                    (sortedAlpha[h] - sortedAlpha[u]);
            final double bScalar = (sortedAlpha[h] - sortedAlpha[h - 1]) /
                    (sortedAlpha[h] - sortedAlpha[u]);

            reverseRecursionRow(
                    curOff + uRow,
                    curOffKp1 + uRow,
                    prevOff,
                    uRow,
                    aScalar,
                    bScalar,
                    transitionMatrix,
                    C,
                    sortedAlpha,
                    h,
                    u,
                    true,
                    collectTransitionAdjoints,
                    collectRewardRateAdjoints);
        }
    }

    private void reverseHighRows(
            int h,
            int n,
            int k,
            double[] transitionMatrix,
            double[] sortedAlpha,
            SericolaCumulantMatrices cumulants,
            boolean collectTransitionAdjoints,
            boolean collectRewardRateAdjoints) {

        final double[] C = cumulants.values();
        final int curOff = cumulants.offset(h, n, k);
        final int curOffKm1 = cumulants.offset(h, n, k - 1);
        final int prevOff = cumulants.offset(h, n - 1, k - 1);

        for (int u = h; u <= phi; ++u) {
            final int uRow = u * dim;
            final double cScalar = (sortedAlpha[u] - sortedAlpha[h]) /
                    (sortedAlpha[u] - sortedAlpha[h - 1]);
            final double dScalar = (sortedAlpha[h] - sortedAlpha[h - 1]) /
                    (sortedAlpha[u] - sortedAlpha[h - 1]);

            reverseRecursionRow(
                    curOff + uRow,
                    curOffKm1 + uRow,
                    prevOff,
                    uRow,
                    cScalar,
                    dScalar,
                    transitionMatrix,
                    C,
                    sortedAlpha,
                    h,
                    u,
                    false,
                    collectTransitionAdjoints,
                    collectRewardRateAdjoints);
        }
    }

    private void reverseRecursionRow(
            int currentRowOffset,
            int adjacentRowOffset,
            int previousOffset,
            int transitionRowOffset,
            double adjacentScale,
            double transitionScale,
            double[] transitionMatrix,
            double[] C,
            double[] sortedAlpha,
            int h,
            int row,
            boolean lowRows,
            boolean collectTransitionAdjoints,
            boolean collectRewardRateAdjoints) {

        boolean hasAdjoint = false;
        for (int v = 0; v < dim; ++v) {
            if (cumulantAdjoints[currentRowOffset + v] != 0.0) {
                hasAdjoint = true;
                break;
            }
        }
        if (!hasAdjoint) {
            return;
        }

        double adjacentScaleAdjoint = 0.0;
        double transitionScaleAdjoint = 0.0;

        for (int w = 0; w < dim; ++w) {
            final int previousRowOffset = previousOffset + w * dim;
            final double transitionEntry = transitionMatrix[transitionRowOffset + w];
            double transitionEntryAdjoint = 0.0;

            for (int v = 0; v < dim; ++v) {
                final double rowAdjoint = cumulantAdjoints[currentRowOffset + v];
                transitionEntryAdjoint += rowAdjoint * C[previousRowOffset + v];
                cumulantAdjoints[previousRowOffset + v] +=
                        transitionScale * transitionEntry * rowAdjoint;
            }

            transitionScaleAdjoint += transitionEntry * transitionEntryAdjoint;
            if (collectTransitionAdjoints) {
                transitionMatrixAdjoint[transitionRowOffset + w] += transitionScale * transitionEntryAdjoint;
            }
        }

        for (int v = 0; v < dim; ++v) {
            final double rowAdjoint = cumulantAdjoints[currentRowOffset + v];
            adjacentScaleAdjoint += rowAdjoint * C[adjacentRowOffset + v];
            cumulantAdjoints[adjacentRowOffset + v] += adjacentScale * rowAdjoint;
            cumulantAdjoints[currentRowOffset + v] = 0.0;
        }

        if (collectRewardRateAdjoints) {
            if (lowRows) {
                accumulateLowRowRewardRateAdjoints(h, row, adjacentScaleAdjoint, transitionScaleAdjoint, sortedAlpha);
            } else {
                accumulateHighRowRewardRateAdjoints(h, row, adjacentScaleAdjoint, transitionScaleAdjoint, sortedAlpha);
            }
        }
    }

    private void accumulateLowRowRewardRateAdjoints(
            int h,
            int row,
            double aAdjoint,
            double bAdjoint,
            double[] sortedAlpha) {

        final double low = sortedAlpha[h - 1];
        final double high = sortedAlpha[h];
        final double rowAlpha = sortedAlpha[row];
        final double width = high - low;
        final double den = high - rowAlpha;
        final double invDen = 1.0 / den;
        final double invDen2 = invDen * invDen;

        rewardRateAdjoint[h - 1] += aAdjoint * invDen - bAdjoint * invDen;
        rewardRateAdjoint[h] +=
                aAdjoint * (-(low - rowAlpha) * invDen2) +
                        bAdjoint * (invDen - width * invDen2);
        rewardRateAdjoint[row] +=
                aAdjoint * (-width * invDen2) +
                        bAdjoint * (width * invDen2);
    }

    private void accumulateHighRowRewardRateAdjoints(
            int h,
            int row,
            double cAdjoint,
            double dAdjoint,
            double[] sortedAlpha) {

        final double low = sortedAlpha[h - 1];
        final double high = sortedAlpha[h];
        final double rowAlpha = sortedAlpha[row];
        final double width = high - low;
        final double den = rowAlpha - low;
        final double invDen = 1.0 / den;
        final double invDen2 = invDen * invDen;

        rewardRateAdjoint[h - 1] +=
                cAdjoint * ((rowAlpha - high) * invDen2) +
                        dAdjoint * (-invDen + width * invDen2);
        rewardRateAdjoint[h] +=
                cAdjoint * (-invDen) +
                        dAdjoint * invDen;
        rewardRateAdjoint[row] +=
                cAdjoint * (width * invDen2) +
                        dAdjoint * (-width * invDen2);
    }

    private void reverseLowBoundaryCopy(int destinationH, int n, SericolaCumulantMatrices cumulants) {
        final int sourceOffset = cumulants.offset(destinationH + 1, n, 0);
        final int destinationOffset = cumulants.offset(destinationH, n, n);

        for (int u = 0; u <= destinationH - 1; ++u) {
            final int row = u * dim;
            for (int v = 0; v < dim; ++v) {
                final int src = sourceOffset + row + v;
                final int dst = destinationOffset + row + v;
                cumulantAdjoints[src] += cumulantAdjoints[dst];
                cumulantAdjoints[dst] = 0.0;
            }
        }
    }

    private void reverseHighBoundaryCopy(int sourceH, int n, SericolaCumulantMatrices cumulants) {
        final int sourceOffset = cumulants.offset(sourceH, n, n);
        final int destinationOffset = cumulants.offset(sourceH + 1, n, 0);

        for (int u = sourceH + 1; u <= phi; ++u) {
            final int row = u * dim;
            for (int v = 0; v < dim; ++v) {
                final int src = sourceOffset + row + v;
                final int dst = destinationOffset + row + v;
                cumulantAdjoints[src] += cumulantAdjoints[dst];
                cumulantAdjoints[dst] = 0.0;
            }
        }
    }

    private void reverseTerminalOverwrite(int n, SericolaCumulantMatrices cumulants, boolean collectTransitionAdjoints) {
        final int terminalOffset = cumulants.offset(phi, n, n);
        final int powerOffset = n * dim2;

        for (int u = 0; u <= phi - 1; ++u) {
            final int row = u * dim;
            for (int v = 0; v < dim; ++v) {
                final int uv = row + v;
                if (collectTransitionAdjoints) {
                    transitionPowerAdjoints[powerOffset + uv] += cumulantAdjoints[terminalOffset + uv];
                }
                cumulantAdjoints[terminalOffset + uv] = 0.0;
            }
        }
    }

    private void computeTransitionPowers(int maxExternal, double[] transitionMatrix) {
        Arrays.fill(transitionPowers, 0, (maxExternal + 1) * dim2, 0.0);
        for (int u = 0; u < dim; ++u) {
            transitionPowers[u * dim + u] = 1.0;
        }

        for (int n = 1; n <= maxExternal; ++n) {
            final int previousOffset = (n - 1) * dim2;
            final int currentOffset = n * dim2;
            rightMultiply(
                    transitionPowers,
                    previousOffset,
                    transitionMatrix,
                    transitionPowers,
                    currentOffset);
        }
    }

    private void reverseTransitionPowers(int maxExternal, double[] transitionMatrix) {
        for (int n = maxExternal; n >= 1; --n) {
            final int previousOffset = (n - 1) * dim2;
            final int currentOffset = n * dim2;

            for (int i = 0; i < dim; ++i) {
                final int iRow = i * dim;
                for (int k = 0; k < dim; ++k) {
                    final int kRow = k * dim;
                    final double previousPower = transitionPowers[previousOffset + iRow + k];
                    for (int j = 0; j < dim; ++j) {
                        final double currentAdjoint = transitionPowerAdjoints[currentOffset + iRow + j];
                        transitionMatrixAdjoint[kRow + j] += previousPower * currentAdjoint;
                        transitionPowerAdjoints[previousOffset + iRow + k] +=
                                currentAdjoint * transitionMatrix[kRow + j];
                    }
                }
            }

            Arrays.fill(transitionPowerAdjoints, currentOffset, currentOffset + dim2, 0.0);
        }
    }

    private void rightMultiply(
            double[] left,
            int leftOffset,
            double[] right,
            double[] out,
            int outOffset) {

        for (int i = 0; i < dim; ++i) {
            final int iRow = i * dim;
            for (int j = 0; j < dim; ++j) {
                double sum = 0.0;
                for (int k = 0; k < dim; ++k) {
                    sum += left[leftOffset + iRow + k] * right[k * dim + j];
                }
                out[outOffset + iRow + j] = sum;
            }
        }
    }

    private void ensureCapacity(int requiredN) {
        if (requiredN < 0) {
            throw new IllegalStateException("C cache must be allocated before gradient evaluation");
        }
        if (requiredN <= allocatedN && cumulantAdjoints != null) {
            return;
        }

        final int n1 = requiredN + 1;
        final long cumulantDoubles = (long) (phi + 1) * (long) n1 * (long) n1 * (long) dim2;
        final long powerDoubles = (long) n1 * (long) dim2;

        if (cumulantDoubles > Integer.MAX_VALUE) {
            throw new IllegalStateException("cumulant adjoint array too large: " + cumulantDoubles + " doubles");
        }
        if (powerDoubles > Integer.MAX_VALUE) {
            throw new IllegalStateException("transition power array too large: " + powerDoubles + " doubles");
        }

        cumulantAdjoints = new double[(int) cumulantDoubles];
        transitionPowers = new double[(int) powerDoubles];
        transitionPowerAdjoints = new double[(int) powerDoubles];
        transitionMatrixAdjoint = new double[dim2];
        sortedDensityAdjoint = new double[dim2];
        rewardRateAdjoint = new double[dim];
        allocatedN = requiredN;
    }

    private void clearAdjoints(int maxExternal) {
        Arrays.fill(cumulantAdjoints, 0, blockEndOffset(maxExternal), 0.0);
        Arrays.fill(transitionMatrixAdjoint, 0.0);
        Arrays.fill(transitionPowerAdjoints, 0, (maxExternal + 1) * dim2, 0.0);
        Arrays.fill(rewardRateAdjoint, 0.0);
    }

    private int blockEndOffset(int maxExternal) {
        return (((phi * (allocatedN + 1) + maxExternal) * (allocatedN + 1) + maxExternal) * dim2) + dim2;
    }

}
