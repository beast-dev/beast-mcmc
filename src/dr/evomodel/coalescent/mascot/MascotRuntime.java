/*
 * MascotRuntime.java
 *
 * Copyright (c) 2002-2026 Alexei Drummond, Andrew Rambaut and Marc Suchard
 */

package dr.evomodel.coalescent.mascot;

import java.util.Arrays;

/** Package-private runtime state and array helpers for the generic MASCOT delegate. */
final class MascotRuntime {

    static final double TIME_TOLERANCE = 1.0e-14;

    private MascotRuntime() {
    }

    static double[] ensure(double[] array, int size) {
        return (array == null || array.length < size) ? new double[size] : array;
    }

    static int[] ensureInt(int[] array, int size) {
        return (array == null || array.length < size) ? new int[size] : array;
    }

    static void addScaledInto(double[] x, double[] dx, double scale, double[] out, int n) {
        for (int i = 0; i < n; i++) {
            out[i] = x[i] + scale * dx[i];
        }
    }

    static void scaleInto(double[] x, double scale, double[] out, int n) {
        for (int i = 0; i < n; i++) {
            out[i] = scale * x[i];
        }
    }

    static void addInPlace(double[] destination, double[] source, int n) {
        for (int i = 0; i < n; i++) {
            destination[i] += source[i];
        }
    }

    static void addScaledInPlace(double[] destination, double[] source, double scale, int n) {
        for (int i = 0; i < n; i++) {
            destination[i] += scale * source[i];
        }
    }

    static final class ActiveState {
        int[] activeIds;
        double[] probabilities;
        int activeCount;
        double logLikelihood;
        int stateCount;
        int[] lineageToActiveIndex;
        int[] lineageGeneration;
        int currentGeneration;

        ActiveState() {
            this.activeIds = new int[1];
            this.probabilities = new double[1];
            this.lineageToActiveIndex = new int[1];
            this.lineageGeneration = new int[1];
            this.currentGeneration = 1;
        }

        void reset(int stateCount, int capacity, int maxLineageId) {
            this.stateCount = stateCount;
            if (activeIds.length < capacity) {
                activeIds = new int[capacity];
            }
            int probabilityCapacity = capacity * stateCount;
            if (probabilities.length < probabilityCapacity) {
                probabilities = new double[probabilityCapacity];
            }
            int lineageCapacity = maxLineageId + 1;
            if (lineageGeneration.length < lineageCapacity) {
                lineageToActiveIndex = new int[lineageCapacity];
                lineageGeneration = new int[lineageCapacity];
            }
            activeCount = 0;
            logLikelihood = 0.0;
            currentGeneration++;
            if (currentGeneration == Integer.MAX_VALUE) {
                Arrays.fill(lineageGeneration, 0);
                currentGeneration = 1;
            }
        }

        void ensureCapacity(int capacity) {
            if (activeIds.length < capacity) {
                activeIds = Arrays.copyOf(activeIds, capacity);
            }
            int probabilityCapacity = capacity * stateCount;
            if (probabilities.length < probabilityCapacity) {
                probabilities = Arrays.copyOf(probabilities, probabilityCapacity);
            }
        }

        boolean isActive(int lineageId) {
            return lineageId >= 0 && lineageId < lineageGeneration.length
                    && lineageGeneration[lineageId] == currentGeneration;
        }

        int activeIndexOf(int lineageId) {
            return isActive(lineageId) ? lineageToActiveIndex[lineageId] : -1;
        }

        void setActiveIndex(int lineageId, int index) {
            lineageGeneration[lineageId] = currentGeneration;
            lineageToActiveIndex[lineageId] = index;
        }

        void clearActiveIndex(int lineageId) {
            if (lineageId >= 0 && lineageId < lineageGeneration.length) {
                lineageGeneration[lineageId] = 0;
            }
        }

        double[] copyProbabilities() {
            return Arrays.copyOf(probabilities, activeCount * stateCount);
        }

        int[] copyActiveIds() {
            return Arrays.copyOf(activeIds, activeCount);
        }
    }

    static final class Workspace {
        double[] integrationState;
        double[] integrationOut;
        double[] activeClockRates;

        double[] k1;
        double[] k2;
        double[] k3;
        double[] k4;
        double[] y2;
        double[] y3;
        double[] y4;

        double[] sums;
        double[] sumsSquares;
        double[] hValues;
        double[] rValues;
        double[] bValues;
        double[] cValues;
        double[] cSums;
        double[] gradQ;
        double[] migrationGram;

        double[] adjointY0;
        double[] adjointK1;
        double[] adjointK2;
        double[] adjointK3;
        double[] adjointK4;
        double[] vjpY;

        double[] reverseCursorA;
        double[] reverseCursorB;
        double[] reverseOperationA;
        double[] reverseOperationB;

        double[] coalP1;
        double[] coalP2;
        double[] coalParent;
        double[] coalescentTimes;
    }

    static final class EpochRates {
        final double[] migrationMatrix;
        final double[] migrationRates;
        final double[] inversePopulation;

        EpochRates(int stateCount) {
            this.migrationMatrix = new double[stateCount * stateCount];
            this.migrationRates = new double[stateCount * (stateCount - 1)];
            this.inversePopulation = new double[stateCount];
        }
    }
}
