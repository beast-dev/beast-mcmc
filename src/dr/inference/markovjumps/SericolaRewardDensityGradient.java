/*
 * SericolaRewardDensityGradient.java
 *
 * Reverse-mode kernel for reward-density gradients with respect to the
 * uniformization transition matrix.
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
                cumulants);

        computeTransitionPowers(maxExternal, transitionMatrix);
        reverseCumulants(maxExternal, transitionMatrix, sortedAlpha, cumulants);
        reverseTransitionPowers(maxExternal, transitionMatrix);

        return directLambdaAdjoint;
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

            previousPremult = premult;
            premult *= lambdaTime / (n + 1.0);
            if (!xIsZero && !xIsOne) {
                w0 *= oneMinus;
            }
        }

        if (conditionalOnZ0) {
            directLambdaAdjoint -= conditionalObjective * time * Math.exp(-lambdaTime) / denominator;
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

    private void reverseCumulants(
            int maxExternal,
            double[] transitionMatrix,
            double[] sortedAlpha,
            SericolaCumulantMatrices cumulants) {

        for (int n = maxExternal; n >= 1; --n) {
            reverseBackwardSweep(n, transitionMatrix, sortedAlpha, cumulants);
            reverseTerminalOverwrite(n, cumulants);
            reverseForwardSweep(n, transitionMatrix, sortedAlpha, cumulants);
        }
    }

    private void reverseBackwardSweep(
            int n,
            double[] transitionMatrix,
            double[] sortedAlpha,
            SericolaCumulantMatrices cumulants) {

        for (int h = 1; h <= phi; ++h) {
            for (int k = 0; k <= n - 1; ++k) {
                reverseLowRows(h, n, k, transitionMatrix, sortedAlpha, cumulants);
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
            SericolaCumulantMatrices cumulants) {

        for (int h = phi; h >= 1; --h) {
            if (h < phi) {
                reverseHighBoundaryCopy(h, n, cumulants);
            }

            for (int k = n; k >= 1; --k) {
                reverseHighRows(h, n, k, transitionMatrix, sortedAlpha, cumulants);
            }
        }
    }

    private void reverseLowRows(
            int h,
            int n,
            int k,
            double[] transitionMatrix,
            double[] sortedAlpha,
            SericolaCumulantMatrices cumulants) {

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
                    C);
        }
    }

    private void reverseHighRows(
            int h,
            int n,
            int k,
            double[] transitionMatrix,
            double[] sortedAlpha,
            SericolaCumulantMatrices cumulants) {

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
                    C);
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
            double[] C) {

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

            transitionMatrixAdjoint[transitionRowOffset + w] += transitionScale * transitionEntryAdjoint;
        }

        for (int v = 0; v < dim; ++v) {
            final double rowAdjoint = cumulantAdjoints[currentRowOffset + v];
            cumulantAdjoints[adjacentRowOffset + v] += adjacentScale * rowAdjoint;
            cumulantAdjoints[currentRowOffset + v] = 0.0;
        }
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

    private void reverseTerminalOverwrite(int n, SericolaCumulantMatrices cumulants) {
        final int terminalOffset = cumulants.offset(phi, n, n);
        final int powerOffset = n * dim2;

        for (int u = 0; u <= phi - 1; ++u) {
            final int row = u * dim;
            for (int v = 0; v < dim; ++v) {
                final int uv = row + v;
                transitionPowerAdjoints[powerOffset + uv] += cumulantAdjoints[terminalOffset + uv];
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
        allocatedN = requiredN;
    }

    private void clearAdjoints(int maxExternal) {
        Arrays.fill(cumulantAdjoints, 0, blockEndOffset(maxExternal), 0.0);
        Arrays.fill(transitionMatrixAdjoint, 0.0);
        Arrays.fill(transitionPowerAdjoints, 0, (maxExternal + 1) * dim2, 0.0);
    }

    private int blockEndOffset(int maxExternal) {
        return (((phi * (allocatedN + 1) + maxExternal) * (allocatedN + 1) + maxExternal) * dim2) + dim2;
    }

}
