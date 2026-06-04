/*
 * SericolaCumulantMatrices.java
 *
 * Flat, reusable storage for Sericola cumulant matrices and their
 * derivative-time second differences.
 */

package dr.inference.markovjumps;

import java.util.Arrays;

final class SericolaCumulantMatrices {

    private final int dim;
    private final int dim2;
    private final int phi;

    private final double[] powerScratch;
    private final double[] multiplyScratch;

    private double[] cumulants;
    private int allocatedN = -1;
    private int allocatedN1 = 0;
    private int computedN = -1;
    private double maxTime = 0.0;

    private double[] secondDifferences;
    private int[] secondDifferenceStamp;
    private int secondDifferenceEpoch = 1;
    private int secondDifferenceAllocatedN = -1;

    SericolaCumulantMatrices(int dim) {
        this.dim = dim;
        this.dim2 = dim * dim;
        this.phi = dim - 1;
        this.powerScratch = new double[dim2];
        this.multiplyScratch = new double[dim2];
    }

    void ensureForTime(double time, int requiredN, double[] transitionMatrix, double[] sortedAlpha) {
        if (time <= 0.0) {
            throw new IllegalArgumentException("time must be > 0");
        }

        ensureCapacityN(requiredN);

        if (time > maxTime) {
            maxTime = time;
        }

        if (computedN < requiredN) {
            computeCumulants(requiredN, transitionMatrix, sortedAlpha);
            computedN = requiredN;
        }
    }

    double[] values() {
        return cumulants;
    }

    int computedN() {
        return computedN;
    }

    int allocatedN() {
        return allocatedN;
    }

    double maxTime() {
        return maxTime;
    }

    int offset(int h, int n, int k) {
        return (((h * allocatedN1 + n) * allocatedN1 + k) * dim2);
    }

    void invalidateComputedExtent() {
        computedN = -1;
        maxTime = 0.0;
        invalidateSecondDifferenceCache();
    }

    void ensureSecondDifferenceCapacity() {
        if (allocatedN < 0) {
            throw new IllegalStateException("C cache must be allocated first");
        }

        if (secondDifferences != null && secondDifferenceAllocatedN == allocatedN) {
            return;
        }

        final long nBlocks = (long) (phi + 1) * (long) allocatedN1 * (long) allocatedN1;
        final long totalDoubles = nBlocks * (long) dim2;

        if (totalDoubles > Integer.MAX_VALUE) {
            throw new IllegalStateException("D2 array too large: " + totalDoubles + " doubles");
        }
        if (nBlocks > Integer.MAX_VALUE) {
            throw new IllegalStateException("d2Stamp array too large: " + nBlocks + " blocks");
        }

        secondDifferences = new double[(int) totalDoubles];
        secondDifferenceStamp = new int[(int) nBlocks];
        secondDifferenceAllocatedN = allocatedN;
        secondDifferenceEpoch = 1;
    }

    void prepareSecondDifferenceRow(int h, int n) {
        final int n1 = n + 1;
        final int baseBlock = secondDifferenceBlockIndex(h, n, 0);
        final int baseOffset = baseBlock * dim2;
        final double[] C = cumulants;
        final double[] D2 = secondDifferences;

        for (int k = 0; k <= n - 1; ++k) {
            final int block = baseBlock + k;
            if (secondDifferenceStamp[block] != secondDifferenceEpoch) {
                final int d2Off = baseOffset + k * dim2;
                final int c0 = offset(h, n1, k);
                final int c1 = offset(h, n1, k + 1);
                final int c2 = offset(h, n1, k + 2);

                for (int uv = 0; uv < dim2; ++uv) {
                    D2[d2Off + uv] = C[c2 + uv] - 2.0 * C[c1 + uv] + C[c0 + uv];
                }
                secondDifferenceStamp[block] = secondDifferenceEpoch;
            }
        }
    }

    int secondDifferenceOffset(int h, int n, int k) {
        return secondDifferenceBlockIndex(h, n, k) * dim2;
    }

    double[] secondDifferences() {
        return secondDifferences;
    }

    private void ensureCapacityN(int requiredN) {
        int newAlloc = (allocatedN < 0) ? requiredN : allocatedN;
        while (newAlloc < requiredN) {
            newAlloc = Math.max(requiredN, (int) (newAlloc * 1.5) + 1);
        }

        if (newAlloc <= allocatedN && cumulants != null) {
            return;
        }

        final int newN1 = newAlloc + 1;
        final long totalDoubles = (long) (phi + 1) * (long) newN1 * (long) newN1 * (long) dim2;

        if (totalDoubles > Integer.MAX_VALUE) {
            throw new IllegalStateException("C array too large: " + totalDoubles + " doubles");
        }

        cumulants = new double[(int) totalDoubles];
        allocatedN = newAlloc;
        allocatedN1 = newN1;

        computedN = -1;
        maxTime = 0.0;
        invalidateSecondDifferenceCache();
    }

    private void computeCumulants(int N, double[] transitionMatrix, double[] sortedAlpha) {
        if (cumulants == null || allocatedN < 0) {
            throw new IllegalStateException("C storage not initialized");
        }

        clearCumulantsUpToN(N);

        Arrays.fill(powerScratch, 0.0);
        for (int u = 0; u < dim; ++u) {
            powerScratch[index(u, u)] = 1.0;
        }

        for (int h = 1; h <= phi; ++h) {
            final int off = offset(h, 0, 0);
            for (int u = 0; u <= h - 1; ++u) {
                cumulants[off + index(u, u)] = 1.0;
            }
        }

        for (int n = 1; n <= N; ++n) {
            computeForwardSweep(n, transitionMatrix, sortedAlpha);

            rightMultiply(powerScratch, transitionMatrix, multiplyScratch);
            System.arraycopy(multiplyScratch, 0, powerScratch, 0, dim2);

            final int off = offset(phi, n, n);
            for (int u = 0; u <= phi - 1; ++u) {
                System.arraycopy(powerScratch, u * dim, cumulants, off + u * dim, dim);
            }

            computeBackwardSweep(n, transitionMatrix, sortedAlpha);
        }
    }

    private void computeForwardSweep(int n, double[] transitionMatrix, double[] sortedAlpha) {
        for (int h = 1; h <= phi; ++h) {
            for (int k = 1; k <= n; ++k) {
                for (int u = h; u <= phi; ++u) {
                    final int uvRow = u * dim;

                    final double cScalar = (sortedAlpha[u] - sortedAlpha[h]) /
                            (sortedAlpha[u] - sortedAlpha[h - 1]);
                    final double dScalar = (sortedAlpha[h] - sortedAlpha[h - 1]) /
                            (sortedAlpha[u] - sortedAlpha[h - 1]);

                    final int curOff = offset(h, n, k);
                    final int curOffKm1 = offset(h, n, k - 1);
                    final int prevOff = offset(h, n - 1, k - 1);

                    for (int v = 0; v <= phi; ++v) {
                        double cVal = cScalar * cumulants[curOffKm1 + (uvRow + v)];

                        double dVal = 0.0;
                        for (int w = 0; w <= phi; ++w) {
                            dVal += transitionMatrix[uvRow + w] * cumulants[prevOff + (w * dim + v)];
                        }
                        dVal *= dScalar;

                        cumulants[curOff + (uvRow + v)] = cVal + dVal;
                    }
                }
            }

            if (h < phi) {
                final int srcOff = offset(h, n, n);
                final int dstOff = offset(h + 1, n, 0);
                for (int u = h + 1; u <= phi; ++u) {
                    System.arraycopy(cumulants, srcOff + u * dim, cumulants, dstOff + u * dim, dim);
                }
            }
        }
    }

    private void computeBackwardSweep(int n, double[] transitionMatrix, double[] sortedAlpha) {
        for (int h = phi; h >= 1; --h) {
            for (int k = n - 1; k >= 0; --k) {
                for (int u = 0; u <= h - 1; u++) {
                    final int uvRow = u * dim;

                    final double cScalar = (sortedAlpha[h - 1] - sortedAlpha[u]) /
                            (sortedAlpha[h] - sortedAlpha[u]);
                    final double dScalar = (sortedAlpha[h] - sortedAlpha[h - 1]) /
                            (sortedAlpha[h] - sortedAlpha[u]);

                    final int curOff = offset(h, n, k);
                    final int curOffKp1 = offset(h, n, k + 1);
                    final int prevOff = offset(h, n - 1, k);

                    for (int v = 0; v <= phi; ++v) {
                        double cVal = cScalar * cumulants[curOffKp1 + (uvRow + v)];

                        double dVal = 0.0;
                        for (int w = 0; w <= phi; ++w) {
                            dVal += transitionMatrix[uvRow + w] * cumulants[prevOff + (w * dim + v)];
                        }
                        dVal *= dScalar;

                        cumulants[curOff + (uvRow + v)] = cVal + dVal;
                    }
                }

                if (h >= 2) {
                    final int srcOff = offset(h, n, 0);
                    final int dstOff = offset(h - 1, n, n);
                    for (int u = 0; u <= h - 2; ++u) {
                        System.arraycopy(cumulants, srcOff + u * dim, cumulants, dstOff + u * dim, dim);
                    }
                }
            }
        }
    }

    private void rightMultiply(double[] A, double[] B, double[] out) {
        Arrays.fill(out, 0.0);
        for (int i = 0; i < dim; ++i) {
            final int ioff = i * dim;
            for (int k = 0; k < dim; ++k) {
                final double aik = A[ioff + k];
                final int koff = k * dim;
                for (int j = 0; j < dim; ++j) {
                    out[ioff + j] += aik * B[koff + j];
                }
            }
        }
    }

    private int secondDifferenceBlockIndex(int h, int n, int k) {
        return ((h * allocatedN1 + n) * allocatedN1 + k);
    }

    private void clearCumulantsUpToN(int N) {
        final int end = offset(phi, N, N) + dim2;
        Arrays.fill(cumulants, 0, end, 0.0);
    }

    private void invalidateSecondDifferenceCache() {
        if (secondDifferenceStamp == null) {
            return;
        }
        secondDifferenceEpoch++;
        if (secondDifferenceEpoch == Integer.MAX_VALUE) {
            Arrays.fill(secondDifferenceStamp, 0);
            secondDifferenceEpoch = 1;
        }
    }

    private int index(int i, int j) {
        return i * dim + j;
    }
}
