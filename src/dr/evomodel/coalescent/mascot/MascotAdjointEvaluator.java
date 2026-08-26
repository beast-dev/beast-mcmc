/*
 * MascotAdjointEvaluator.java
 *
 * Copyright (c) 2002-2026 Alexei Drummond, Andrew Rambaut and Marc Suchard
 */

package dr.evomodel.coalescent.mascot;

import java.util.Arrays;

/** Owns MASCOT adjoint recording storage and reverse-mode replay. */
final class MascotAdjointEvaluator {

    private static final double EPS = 1.0e-300;

    private final int stateCount;
    private final int parametersPerEpoch;
    private final int migrationParametersPerEpoch;
    private final int parameterCount;
    private final MascotAdjointTape.Store operationTape = new MascotAdjointTape.Store();
    private double[] gradientScratch;
    private double[] clockGradientScratch;

    MascotAdjointEvaluator(int stateCount, int parametersPerEpoch, int migrationParametersPerEpoch,
                            int parameterCount) {
        this.stateCount = stateCount;
        this.parametersPerEpoch = parametersPerEpoch;
        this.migrationParametersPerEpoch = migrationParametersPerEpoch;
        this.parameterCount = parameterCount;
    }

    MascotAdjointTape.Store resetTape(int expectedOperationCount) {
        operationTape.reset(expectedOperationCount);
        return operationTape;
    }

    void reverseInto(MascotRuntime.Workspace workspace, MascotRuntime.EpochRates[] epochRates,
                     MascotAdjointTape.Store operations, int finalActiveCount,
                     int lineageDimension, boolean needsClockGradientScratch,
                     double[] gradientOut, double[] clockGradientOut, double[] ancestralStateScores) {
        double[] gradient = gradientOut;
        if (gradient == null) {
            gradientScratch = MascotRuntime.ensure(gradientScratch, parameterCount);
            gradient = gradientScratch;
        }
        double[] clockGradient = clockGradientOut;
        if (clockGradient == null && needsClockGradientScratch) {
            clockGradientScratch = MascotRuntime.ensure(clockGradientScratch, lineageDimension);
            clockGradient = clockGradientScratch;
        }
        reverse(workspace, epochRates, operations, finalActiveCount, gradient, clockGradient, ancestralStateScores);
    }

    private void reverse(MascotRuntime.Workspace workspace, MascotRuntime.EpochRates[] epochRates,
                         MascotAdjointTape.Store operations, int finalActiveCount, double[] gradient,
                         double[] clockGradient, double[] ancestralStateScores) {
        Arrays.fill(gradient, 0, parameterCount, 0.0);
        if (clockGradient != null) {
            Arrays.fill(clockGradient, 0, clockGradient.length, 0.0);
        }
        int dim = finalActiveCount * stateCount + 1;

        workspace.reverseOperationA = MascotRuntime.ensure(workspace.reverseOperationA, dim);
        double[] cursor = workspace.reverseOperationA;
        Arrays.fill(cursor, 0, dim, 0.0);
        cursor[dim - 1] = 1.0;
        boolean cursorIsA = true;

        for (int opIndex = operations.size() - 1; opIndex >= 0; opIndex--) {
            MascotAdjointTape.Operation operation = operations.get(opIndex);
            int nextDim;
            if (operation instanceof MascotAdjointTape.Interval) {
                nextDim = dim;
            } else if (operation instanceof MascotAdjointTape.Coalescent) {
                nextDim = dim + stateCount;
            } else if (operation instanceof MascotAdjointTape.Sample) {
                nextDim = dim - stateCount;
            } else {
                throw new IllegalArgumentException("unknown tape operation: " + operation.getClass());
            }

            if (operation instanceof MascotAdjointTape.Sample) {
                reverseSampleInPlace((MascotAdjointTape.Sample) operation, cursor, dim);
                dim = nextDim;
                continue;
            }

            double[] next;
            if (cursorIsA) {
                workspace.reverseOperationB = MascotRuntime.ensure(workspace.reverseOperationB, nextDim);
                next = workspace.reverseOperationB;
            } else {
                workspace.reverseOperationA = MascotRuntime.ensure(workspace.reverseOperationA, nextDim);
                next = workspace.reverseOperationA;
            }

            if (operation instanceof MascotAdjointTape.Interval) {
                reverseIntervalInto(workspace, (MascotAdjointTape.Interval) operation, cursor, next, epochRates,
                        gradient, clockGradient);
            } else if (operation instanceof MascotAdjointTape.Coalescent) {
                reverseCoalescentInto((MascotAdjointTape.Coalescent) operation, cursor, dim, next, gradient,
                        ancestralStateScores);
            }

            cursor = next;
            cursorIsA = !cursorIsA;
            dim = nextDim;
        }
    }

    private void reverseIntervalInto(MascotRuntime.Workspace workspace, MascotAdjointTape.Interval tape,
                                     double[] adjointAfter, double[] adjointBeforeOut,
                                     MascotRuntime.EpochRates[] epochRates, double[] gradient, double[] clockGradient) {
        int dim = tape.stateDimension;
        MascotRuntime.EpochRates rates = epochRates[tape.epoch];

        workspace.reverseCursorA = MascotRuntime.ensure(workspace.reverseCursorA, dim);
        workspace.reverseCursorB = MascotRuntime.ensure(workspace.reverseCursorB, dim);
        System.arraycopy(adjointAfter, 0, workspace.reverseCursorA, 0, dim);
        double[] cursor = workspace.reverseCursorA;
        double[] next = workspace.reverseCursorB;

        for (int step = tape.steps - 1; step >= 0; step--) {
            int offset = step * dim;
            reverseStepInto(tape, offset, rates, cursor, next, gradient, clockGradient, workspace);
            double[] swap = cursor;
            cursor = next;
            next = swap;
        }

        System.arraycopy(cursor, 0, adjointBeforeOut, 0, dim);
    }

    private void reverseStepInto(MascotAdjointTape.Interval tape, int offset, MascotRuntime.EpochRates rates,
                                 double[] adjointAfter, double[] adjointBeforeOut, double[] gradient,
                                 double[] clockGradient, MascotRuntime.Workspace w) {
        int dim = tape.stateDimension;
        int activeCount = tape.activeCount;
        int epoch = tape.epoch;
        double h = tape.h;
        int[] activeIds = tape.activeIds;
        double[] activeClockRates = tape.clockRates;

        w.adjointY0 = MascotRuntime.ensure(w.adjointY0, dim);
        System.arraycopy(adjointAfter, 0, w.adjointY0, 0, dim);

        w.adjointK1 = MascotRuntime.ensure(w.adjointK1, dim);
        MascotRuntime.scaleInto(adjointAfter, h / 6.0, w.adjointK1, dim);
        w.adjointK2 = MascotRuntime.ensure(w.adjointK2, dim);
        MascotRuntime.scaleInto(adjointAfter, h / 3.0, w.adjointK2, dim);
        w.adjointK3 = MascotRuntime.ensure(w.adjointK3, dim);
        MascotRuntime.scaleInto(adjointAfter, h / 3.0, w.adjointK3, dim);
        w.adjointK4 = MascotRuntime.ensure(w.adjointK4, dim);
        MascotRuntime.scaleInto(adjointAfter, h / 6.0, w.adjointK4, dim);

        w.vjpY = MascotRuntime.ensure(w.vjpY, dim);

        rhsVjpInto(tape.y4, offset, activeCount, rates, epoch, w.adjointK4, w.vjpY, gradient, w,
                activeIds, activeClockRates, clockGradient);
        MascotRuntime.addInPlace(w.adjointY0, w.vjpY, dim);
        MascotRuntime.addScaledInPlace(w.adjointK3, w.vjpY, h, dim);

        rhsVjpInto(tape.y3, offset, activeCount, rates, epoch, w.adjointK3, w.vjpY, gradient, w,
                activeIds, activeClockRates, clockGradient);
        MascotRuntime.addInPlace(w.adjointY0, w.vjpY, dim);
        MascotRuntime.addScaledInPlace(w.adjointK2, w.vjpY, 0.5 * h, dim);

        rhsVjpInto(tape.y2, offset, activeCount, rates, epoch, w.adjointK2, w.vjpY, gradient, w,
                activeIds, activeClockRates, clockGradient);
        MascotRuntime.addInPlace(w.adjointY0, w.vjpY, dim);
        MascotRuntime.addScaledInPlace(w.adjointK1, w.vjpY, 0.5 * h, dim);

        rhsVjpInto(tape.y0, offset, activeCount, rates, epoch, w.adjointK1, w.vjpY, gradient, w,
                activeIds, activeClockRates, clockGradient);
        MascotRuntime.addInPlace(w.adjointY0, w.vjpY, dim);

        System.arraycopy(w.adjointY0, 0, adjointBeforeOut, 0, dim);
    }

    private void rhsVjpInto(double[] y, int yOffset, int activeCount, MascotRuntime.EpochRates rates, int epoch,
                            double[] adjointRhs, double[] adjointYOut, double[] gradient, MascotRuntime.Workspace w,
                            int[] activeLineageIds, double[] activeClockRates, double[] clockGradient) {
        int K = stateCount;
        int stateSize = activeCount * K;
        adjointYOut[stateSize] = 0.0;

        w.migrationGram = MascotRuntime.ensure(w.migrationGram, K * K);
        w.sums = MascotRuntime.ensure(w.sums, K);
        w.sumsSquares = MascotRuntime.ensure(w.sumsSquares, K);

        {
            double c = activeClockRates == null ? 1.0 : activeClockRates[0];
            double clockContribution = 0.0;
            for (int source = 0; source < K; source++) {
                double ySource = y[yOffset + source];
                w.sums[source] = ySource;
                w.sumsSquares[source] = ySource * ySource;
                double v = 0.0;
                int row = source * K;
                for (int sink = 0; sink < K; sink++) {
                    double adjSink = adjointRhs[sink];
                    v += adjSink * rates.migrationMatrix[row + sink];
                    w.migrationGram[row + sink] = c * ySource * adjSink;
                }
                adjointYOut[source] = c * v;
                if (clockGradient != null) {
                    clockContribution += ySource * v;
                }
            }
            if (clockGradient != null) {
                clockGradient[activeLineageIds[0]] += clockContribution;
            }
        }

        for (int lineage = 1; lineage < activeCount; lineage++) {
            int offset = lineage * K;
            double c = activeClockRates == null ? 1.0 : activeClockRates[lineage];
            double clockContribution = 0.0;
            for (int source = 0; source < K; source++) {
                double ySource = y[yOffset + offset + source];
                w.sums[source] += ySource;
                w.sumsSquares[source] += ySource * ySource;
                double v = 0.0;
                int row = source * K;
                for (int sink = 0; sink < K; sink++) {
                    double adjSink = adjointRhs[offset + sink];
                    v += adjSink * rates.migrationMatrix[row + sink];
                    w.migrationGram[row + sink] += c * ySource * adjSink;
                }
                adjointYOut[offset + source] = c * v;
                if (clockGradient != null) {
                    clockContribution += ySource * v;
                }
            }
            if (clockGradient != null) {
                clockGradient[activeLineageIds[lineage]] += clockContribution;
            }
        }

        int thetaOffset = epoch * parametersPerEpoch;
        int rateIndex = 0;
        for (int source = 0; source < K; source++) {
            int row = source * K;
            double diagonal = w.migrationGram[row + source];
            for (int sink = 0; sink < K; sink++) {
                if (source == sink) {
                    continue;
                }
                double contribution = w.migrationGram[row + sink] - diagonal;
                gradient[thetaOffset + rateIndex] += contribution;
                rateIndex++;
            }
        }

        w.hValues = MascotRuntime.ensure(w.hValues, stateSize);
        w.rValues = MascotRuntime.ensure(w.rValues, activeCount);
        for (int lineage = 0; lineage < activeCount; lineage++) {
            int offset = lineage * K;
            double r = 0.0;
            for (int state = 0; state < K; state++) {
                double h = (w.sums[state] - y[yOffset + offset + state]) * rates.inversePopulation[state];
                w.hValues[offset + state] = h;
                r += y[yOffset + offset + state] * h;
            }
            w.rValues[lineage] = r;
        }

        w.bValues = MascotRuntime.ensure(w.bValues, activeCount);
        for (int lineage = 0; lineage < activeCount; lineage++) {
            int offset = lineage * K;
            double b = 0.0;
            for (int state = 0; state < K; state++) {
                b += adjointRhs[offset + state] * y[yOffset + offset + state];
            }
            w.bValues[lineage] = b;
        }

        w.cSums = MascotRuntime.ensure(w.cSums, K);
        w.cValues = MascotRuntime.ensure(w.cValues, stateSize);
        double b = w.bValues[0];
        for (int state = 0; state < K; state++) {
            double upstream = adjointRhs[state];
            adjointYOut[state] += upstream * (w.rValues[0] - w.hValues[state]) +
                    b * w.hValues[state];
            double c = y[yOffset + state] * (b - upstream);
            w.cValues[state] = c;
            w.cSums[state] = c;
        }
        for (int lineage = 1; lineage < activeCount; lineage++) {
            int offset = lineage * K;
            b = w.bValues[lineage];
            for (int state = 0; state < K; state++) {
                double upstream = adjointRhs[offset + state];
                adjointYOut[offset + state] += upstream * (w.rValues[lineage] - w.hValues[offset + state]) +
                        b * w.hValues[offset + state];
                double c = y[yOffset + offset + state] * (b - upstream);
                w.cValues[offset + state] = c;
                w.cSums[state] += c;
            }
        }

        w.gradQ = MascotRuntime.ensure(w.gradQ, K);
        double ellAdjoint = adjointRhs[stateSize];
        for (int state = 0; state < K; state++) {
            double pairSum = 0.5 * (w.sums[state] * w.sums[state] - w.sumsSquares[state]);
            w.gradQ[state] = -ellAdjoint * pairSum;
        }
        for (int lineage = 0; lineage < activeCount; lineage++) {
            int offset = lineage * K;
            for (int state = 0; state < K; state++) {
                double c = w.cValues[offset + state];
                adjointYOut[offset + state] += rates.inversePopulation[state] * (w.cSums[state] - c) -
                        ellAdjoint * w.hValues[offset + state];
                w.gradQ[state] += c * (w.sums[state] - y[yOffset + offset + state]);
            }
        }

        int etaOffset = thetaOffset + migrationParametersPerEpoch;
        for (int state = 0; state < K; state++) {
            gradient[etaOffset + state] += -rates.inversePopulation[state] * w.gradQ[state];
        }
    }

    private void reverseSampleInPlace(MascotAdjointTape.Sample tape, double[] adjointAfter, int afterDim) {
        int afterCount = (afterDim - 1) / stateCount;
        int beforeCount = afterCount - 1;
        if (tape.sampleIndexAfter != beforeCount) {
            throw new IllegalStateException("sample was not appended at the final active slot");
        }
        int beforeDim = beforeCount * stateCount + 1;
        adjointAfter[beforeDim - 1] = adjointAfter[afterDim - 1];
    }

    private void reverseCoalescentInto(MascotAdjointTape.Coalescent tape, double[] adjointAfter, int afterDim,
                                       double[] adjointBeforeOut, double[] gradient,
                                       double[] ancestralStateScores) {
        int beforeDim = afterDim + stateCount;
        int beforeCount = (beforeDim - 1) / stateCount;

        for (int beforeIndex = 0; beforeIndex < beforeCount; beforeIndex++) {
            if (beforeIndex == tape.child1Index || beforeIndex == tape.child2Index) {
                continue;
            }
            int afterIndex = beforeIndex == tape.movedFromIndexBefore ? tape.movedToIndexAfter : beforeIndex;
            System.arraycopy(adjointAfter, afterIndex * stateCount,
                    adjointBeforeOut, beforeIndex * stateCount, stateCount);
        }

        int parentOffset = tape.parentIndexAfter * stateCount;
        double dot = 0.0;
        for (int s = 0; s < stateCount; s++) {
            dot += adjointAfter[parentOffset + s] * tape.parentProbabilities[s];
        }

        int thetaOffset = tape.epoch * parametersPerEpoch + migrationParametersPerEpoch;
        double ellAdjoint = adjointAfter[afterDim - 1];
        for (int s = 0; s < stateCount; s++) {
            double p1 = Math.max(tape.p1[s], EPS);
            double p2 = Math.max(tape.p2[s], EPS);
            double centered = adjointAfter[parentOffset + s] - dot;

            double nodeStateScore = tape.parentProbabilities[s] * (centered + ellAdjoint);

            adjointBeforeOut[tape.child1Index * stateCount + s] = nodeStateScore / p1;
            adjointBeforeOut[tape.child2Index * stateCount + s] = nodeStateScore / p2;
            gradient[thetaOffset + s] -= nodeStateScore;

            if (ancestralStateScores != null) {
                ancestralStateScores[tape.parentLineageId * stateCount + s] = nodeStateScore;
            }
        }

        adjointBeforeOut[beforeDim - 1] = ellAdjoint;
    }
}
