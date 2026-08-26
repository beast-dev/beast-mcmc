/*
 * MascotAdjointTape.java
 *
 * Copyright (c) 2002-2026 Alexei Drummond, Andrew Rambaut and Marc Suchard
 */

package dr.evomodel.coalescent.mascot;

import java.util.Arrays;

/** Reusable forward records consumed by the MASCOT adjoint evaluator. */
final class MascotAdjointTape {

    private MascotAdjointTape() {
    }

    interface Operation {
    }

    static final class Store {
        private Operation[] operations = new Operation[16];
        private int operationCount;

        void reset(int expectedOperationCount) {
            ensureOperationCapacity(expectedOperationCount);
            operationCount = 0;
        }

        int size() {
            return operationCount;
        }

        Operation get(int index) {
            return operations[index];
        }

        Interval addInterval(int steps, int activeCount, int stateDimension, int epoch, double h,
                             int[] activeIds, double[] activeClockRates) {
            int index = nextIndex();
            Operation operation = operations[index];
            Interval tape;
            if (operation instanceof Interval) {
                tape = (Interval) operation;
            } else {
                tape = new Interval();
                operations[index] = tape;
            }
            tape.reset(steps, activeCount, stateDimension, epoch, h, activeIds, activeClockRates);
            return tape;
        }

        void addSample(int sampleIndexAfter) {
            int index = nextIndex();
            Operation operation = operations[index];
            Sample tape;
            if (operation instanceof Sample) {
                tape = (Sample) operation;
            } else {
                tape = new Sample();
                operations[index] = tape;
            }
            tape.reset(sampleIndexAfter);
        }

        void addCoalescent(int epoch, int child1Index, int child2Index, int parentIndexAfter,
                           int movedFromIndexBefore, int movedToIndexAfter, int parentLineageId,
                           double[] p1, double[] p2, double[] parentProbabilities,
                           double lambda, int stateCount) {
            int index = nextIndex();
            Operation operation = operations[index];
            Coalescent tape;
            if (operation instanceof Coalescent) {
                tape = (Coalescent) operation;
            } else {
                tape = new Coalescent();
                operations[index] = tape;
            }
            tape.reset(epoch, child1Index, child2Index, parentIndexAfter,
                    movedFromIndexBefore, movedToIndexAfter, parentLineageId,
                    p1, p2, parentProbabilities, lambda, stateCount);
        }

        private int nextIndex() {
            ensureOperationCapacity(operationCount + 1);
            return operationCount++;
        }

        private void ensureOperationCapacity(int capacity) {
            if (operations.length < capacity) {
                int newCapacity = operations.length;
                while (newCapacity < capacity) {
                    newCapacity *= 2;
                }
                operations = Arrays.copyOf(operations, newCapacity);
            }
        }
    }

    static final class Interval implements Operation {
        int steps;
        int activeCount;
        int stateDimension;
        int epoch;
        double h;
        double[] y0;
        double[] y2;
        double[] y3;
        double[] y4;
        int[] activeIds;
        double[] clockRates;

        private void reset(int steps, int activeCount, int stateDimension, int epoch, double h,
                           int[] sourceActiveIds, double[] sourceClockRates) {
            this.steps = steps;
            this.activeCount = activeCount;
            this.stateDimension = stateDimension;
            this.epoch = epoch;
            this.h = h;
            int storageSize = steps * stateDimension;
            y0 = MascotRuntime.ensure(y0, storageSize);
            y2 = MascotRuntime.ensure(y2, storageSize);
            y3 = MascotRuntime.ensure(y3, storageSize);
            y4 = MascotRuntime.ensure(y4, storageSize);
            if (sourceActiveIds != null) {
                activeIds = MascotRuntime.ensureInt(activeIds, activeCount);
                System.arraycopy(sourceActiveIds, 0, activeIds, 0, activeCount);
            } else {
                activeIds = null;
            }
            if (sourceClockRates != null) {
                clockRates = MascotRuntime.ensure(clockRates, activeCount);
                System.arraycopy(sourceClockRates, 0, clockRates, 0, activeCount);
            } else {
                clockRates = null;
            }
        }
    }

    static final class Sample implements Operation {
        int sampleIndexAfter;

        private void reset(int sampleIndexAfter) {
            this.sampleIndexAfter = sampleIndexAfter;
        }
    }

    static final class Coalescent implements Operation {
        int epoch;
        int child1Index;
        int child2Index;
        int parentIndexAfter;
        int movedFromIndexBefore;
        int movedToIndexAfter;
        int parentLineageId;
        double[] p1;
        double[] p2;
        double[] parentProbabilities;
        @SuppressWarnings("unused")
        private double lambda;

        private void reset(int epoch, int child1Index, int child2Index, int parentIndexAfter,
                           int movedFromIndexBefore, int movedToIndexAfter, int parentLineageId,
                           double[] p1, double[] p2, double[] parentProbabilities,
                           double lambda, int stateCount) {
            this.epoch = epoch;
            this.child1Index = child1Index;
            this.child2Index = child2Index;
            this.parentIndexAfter = parentIndexAfter;
            this.movedFromIndexBefore = movedFromIndexBefore;
            this.movedToIndexAfter = movedToIndexAfter;
            this.parentLineageId = parentLineageId;
            this.lambda = lambda;
            this.p1 = MascotRuntime.ensure(this.p1, stateCount);
            this.p2 = MascotRuntime.ensure(this.p2, stateCount);
            this.parentProbabilities = MascotRuntime.ensure(this.parentProbabilities, stateCount);
            System.arraycopy(p1, 0, this.p1, 0, stateCount);
            System.arraycopy(p2, 0, this.p2, 0, stateCount);
            System.arraycopy(parentProbabilities, 0, this.parentProbabilities, 0, stateCount);
        }
    }
}
