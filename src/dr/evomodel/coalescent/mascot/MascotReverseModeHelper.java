/*
 * MascotReverseModeHelper.java
 *
 * Copyright (c) 2002-2026 Alexei Drummond, Andrew Rambaut and Marc Suchard
 */

package dr.evomodel.coalescent.mascot;

import java.util.Arrays;

/**
 * The backward (adjoint/reverse-mode) replay over one {@link MascotCore}
 * instance's recorded {@link MascotCore.OperationTapeStore}. {@link
 * MascotCore} itself owns the forward RK4 integration (which also records the
 * tape this class replays -- taping is inherently a forward-pass concern,
 * since the intermediate values it captures can't be reconstructed later) and
 * the {@code evaluate(...)}/{@code evaluateInto(...)} orchestration; this
 * class owns only the backward traversal.
 * <p/>
 * Extracted purely for readability, not as an attempt at a clean forward/reverse
 * architectural boundary: the two passes remain tightly coupled by design (the
 * shared {@link MascotCore.EpochRates} cache, the tape itself, and {@link
 * MascotCore.Workspace}'s reused scratch arrays all cross this class boundary
 * on every call). Every per-evaluation input this class needs -- the {@link
 * MascotCore.Workspace}, the {@link MascotCore.EpochRates} array, and the
 * {@link MascotCore.OperationTapeStore} -- is therefore passed explicitly into
 * {@link #reverse}, rather than this class reaching back into a shared
 * {@code MascotCore} instance for them.
 * <p/>
 * Constructed once per {@code MascotCore} instance with that instance's fixed
 * dimensional constants (mirroring {@code MascotCore}'s own fields). Not
 * thread-safe to call concurrently against the same {@code Workspace}, for the
 * same reason {@link MascotCore} itself isn't (see its class doc).
 */
final class MascotReverseModeHelper {

    private static final double EPS = 1.0e-300;

    private final int stateCount;
    private final int parametersPerEpoch;
    private final int migrationParametersPerEpoch;
    private final int parameterCount;

    MascotReverseModeHelper(int stateCount, int parametersPerEpoch, int migrationParametersPerEpoch,
                            int parameterCount) {
        this.stateCount = stateCount;
        this.parametersPerEpoch = parametersPerEpoch;
        this.migrationParametersPerEpoch = migrationParametersPerEpoch;
        this.parameterCount = parameterCount;
    }

    /**
     * Writes into {@code gradient} (caller/scratch-owned, length {@code
     * parameterCount}) instead of allocating one: fixed-tree HMC calls this
     * once per gradient evaluation, so a fresh array here would be the single
     * largest per-call allocation in the reverse pass.
     */
    void reverse(MascotCore.Workspace workspace, MascotCore.EpochRates[] epochRates,
                MascotCore.OperationTapeStore operations, int finalActiveCount, double[] gradient,
                double[] clockGradient, double[] ancestralStateScores) {
        Arrays.fill(gradient, 0, parameterCount, 0.0);
        if (clockGradient != null) {
            Arrays.fill(clockGradient, 0, clockGradient.length, 0.0);
        }
        int dim = finalActiveCount * stateCount + 1;

        // Ping-pong between two reusable buffers across the whole reverse traversal
        // (operations included) instead of allocating a fresh adjoint array at every
        // sample/coalescent/interval boundary. Dimension changes by +/- stateCount at
        // sample/coalescent events, so each buffer is grown (never shrunk) to the
        // largest dimension it is ever asked to hold. These are distinct from
        // Workspace.reverseCursorA/B, which reverseIntervalInto uses internally for
        // its own per-RK4-step ping-pong (same dimension throughout one interval, so
        // that inner loop does not need to touch these outer, operation-level buffers).
        workspace.reverseOperationA = MascotCore.ensure(workspace.reverseOperationA, dim);
        double[] cursor = workspace.reverseOperationA;
        Arrays.fill(cursor, 0, dim, 0.0);
        cursor[dim - 1] = 1.0;
        boolean cursorIsA = true;

        for (int opIndex = operations.size() - 1; opIndex >= 0; opIndex--) {
            MascotCore.OperationTape operation = operations.get(opIndex);
            int nextDim;
            if (operation instanceof MascotCore.IntervalTape) {
                nextDim = dim;
            } else if (operation instanceof MascotCore.CoalescentTape) {
                nextDim = dim + stateCount;
            } else if (operation instanceof MascotCore.SampleTape) {
                nextDim = dim - stateCount;
            } else {
                throw new IllegalArgumentException("unknown tape operation: " + operation.getClass());
            }

            if (operation instanceof MascotCore.SampleTape) {
                reverseSampleInPlace((MascotCore.SampleTape) operation, cursor, dim);
                dim = nextDim;
                continue;
            }

            double[] next;
            if (cursorIsA) {
                workspace.reverseOperationB = MascotCore.ensure(workspace.reverseOperationB, nextDim);
                next = workspace.reverseOperationB;
            } else {
                workspace.reverseOperationA = MascotCore.ensure(workspace.reverseOperationA, nextDim);
                next = workspace.reverseOperationA;
            }

            if (operation instanceof MascotCore.IntervalTape) {
                reverseIntervalInto(workspace, (MascotCore.IntervalTape) operation, cursor, next, epochRates,
                        gradient, clockGradient);
            } else if (operation instanceof MascotCore.CoalescentTape) {
                reverseCoalescentInto((MascotCore.CoalescentTape) operation, cursor, dim, next, gradient,
                        ancestralStateScores);
            }

            cursor = next;
            cursorIsA = !cursorIsA;
            dim = nextDim;
        }
    }

    private void reverseIntervalInto(MascotCore.Workspace workspace, MascotCore.IntervalTape tape,
                                     double[] adjointAfter, double[] adjointBeforeOut,
                                     MascotCore.EpochRates[] epochRates, double[] gradient, double[] clockGradient) {
        int dim = tape.stateDimension;
        MascotCore.EpochRates rates = epochRates[tape.epoch];

        workspace.reverseCursorA = MascotCore.ensure(workspace.reverseCursorA, dim);
        workspace.reverseCursorB = MascotCore.ensure(workspace.reverseCursorB, dim);
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

    private void reverseStepInto(MascotCore.IntervalTape tape, int offset, MascotCore.EpochRates rates,
                                 double[] adjointAfter, double[] adjointBeforeOut, double[] gradient,
                                 double[] clockGradient, MascotCore.Workspace w) {
        int dim = tape.stateDimension;
        int activeCount = tape.activeCount;
        int epoch = tape.epoch;
        double h = tape.h;
        int[] activeIds = tape.activeIds;
        double[] activeClockRates = tape.clockRates;

        w.adjointY0 = MascotCore.ensure(w.adjointY0, dim);
        System.arraycopy(adjointAfter, 0, w.adjointY0, 0, dim);

        w.adjointK1 = MascotCore.ensure(w.adjointK1, dim);
        MascotCore.scaleInto(adjointAfter, h / 6.0, w.adjointK1, dim);
        w.adjointK2 = MascotCore.ensure(w.adjointK2, dim);
        MascotCore.scaleInto(adjointAfter, h / 3.0, w.adjointK2, dim);
        w.adjointK3 = MascotCore.ensure(w.adjointK3, dim);
        MascotCore.scaleInto(adjointAfter, h / 3.0, w.adjointK3, dim);
        w.adjointK4 = MascotCore.ensure(w.adjointK4, dim);
        MascotCore.scaleInto(adjointAfter, h / 6.0, w.adjointK4, dim);

        w.vjpY = MascotCore.ensure(w.vjpY, dim);

        rhsVjpInto(tape.y4, offset, activeCount, rates, epoch, w.adjointK4, w.vjpY, gradient, w,
                activeIds, activeClockRates, clockGradient);
        MascotCore.addInPlace(w.adjointY0, w.vjpY, dim);
        MascotCore.addScaledInPlace(w.adjointK3, w.vjpY, h, dim);

        rhsVjpInto(tape.y3, offset, activeCount, rates, epoch, w.adjointK3, w.vjpY, gradient, w,
                activeIds, activeClockRates, clockGradient);
        MascotCore.addInPlace(w.adjointY0, w.vjpY, dim);
        MascotCore.addScaledInPlace(w.adjointK2, w.vjpY, 0.5 * h, dim);

        rhsVjpInto(tape.y2, offset, activeCount, rates, epoch, w.adjointK2, w.vjpY, gradient, w,
                activeIds, activeClockRates, clockGradient);
        MascotCore.addInPlace(w.adjointY0, w.vjpY, dim);
        MascotCore.addScaledInPlace(w.adjointK1, w.vjpY, 0.5 * h, dim);

        rhsVjpInto(tape.y0, offset, activeCount, rates, epoch, w.adjointK1, w.vjpY, gradient, w,
                activeIds, activeClockRates, clockGradient);
        MascotCore.addInPlace(w.adjointY0, w.vjpY, dim);

        System.arraycopy(w.adjointY0, 0, adjointBeforeOut, 0, dim);
    }

    /**
     * {@code y} is read at {@code yOffset + ...} throughout, so callers can pass a
     * taped stage array (e.g. {@code tape.y4}) directly at the right per-step
     * offset instead of copying that step's slice into a scratch buffer first.
     * {@code adjointRhs}, {@code adjointYOut}, and {@code gradient} are always
     * 0-based (they are per-call workspace/output buffers, not taped arrays).
     */
    private void rhsVjpInto(double[] y, int yOffset, int activeCount, MascotCore.EpochRates rates, int epoch,
                            double[] adjointRhs, double[] adjointYOut, double[] gradient, MascotCore.Workspace w,
                            int[] activeLineageIds, double[] activeClockRates, double[] clockGradient) {
        int K = stateCount;
        int stateSize = activeCount * K;
        // adjointYOut[stateSize] (the VJP component for the input's log-likelihood
        // slot) is always zero: the forward RHS's out[stateSize] = -hazard never
        // reads y[stateSize], so nothing below ever writes this index again.
        adjointYOut[stateSize] = 0.0;

        // migrationGram[source*K+sink] = sum_lineage clock[lineage] * y[.,source] * adjointRhs[.,sink]
        // (unweighted, i.e. clock[lineage]=1, when clockGradient/activeClockRates are
        // null). The theta-gradient contribution for migration rate (source -> sink) is
        // sum_lineage clock[lineage] * y[.,source] * (adjointRhs[.,sink] - adjointRhs[.,source]),
        // which is exactly migrationGram[source,sink] - migrationGram[source,source]: the
        // diagonal entries are precisely the term the old two-separate-loops version
        // recomputed redundantly (K-1 times per source, once per sink). Building this
        // matrix inside the same lineage-outer pass as the y-adjoint below removes
        // that redundancy and turns the old (source,sink)-outer / lineage-inner
        // stride-K loop into one extra accumulation per already-visited element.
        w.migrationGram = MascotCore.ensure(w.migrationGram, K * K);
        w.sums = MascotCore.ensure(w.sums, K);
        w.sumsSquares = MascotCore.ensure(w.sumsSquares, K);

        // sums/sumsSquares accumulation is fused in here too (same lineage-outer
        // pass, reusing ySource). The first lineage assigns every scratch slot
        // that used to be zero-filled; remaining lineages accumulate into it.
        // sums/sumsSquares are never clock-weighted (see rhsInto's matching note).
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
                    // d(logL)/d(clock[lineage]) accumulates the *unscaled* v here
                    // (the migration RHS value the clock would have multiplied),
                    // per the chain rule derived for out[l,sink] = c_l * (Q^T y)[sink].
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

        w.hValues = MascotCore.ensure(w.hValues, stateSize);
        w.rValues = MascotCore.ensure(w.rValues, activeCount);
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

        w.bValues = MascotCore.ensure(w.bValues, activeCount);
        for (int lineage = 0; lineage < activeCount; lineage++) {
            int offset = lineage * K;
            double b = 0.0;
            for (int state = 0; state < K; state++) {
                b += adjointRhs[offset + state] * y[yOffset + offset + state];
            }
            w.bValues[lineage] = b;
        }

        w.cSums = MascotCore.ensure(w.cSums, K);
        w.cValues = MascotCore.ensure(w.cValues, stateSize);
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

        w.gradQ = MascotCore.ensure(w.gradQ, K);
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

    private void reverseSampleInPlace(MascotCore.SampleTape tape, double[] adjointAfter, int afterDim) {
        int afterCount = (afterDim - 1) / stateCount;
        int beforeCount = afterCount - 1;
        if (tape.sampleIndexAfter != beforeCount) {
            throw new IllegalStateException("sample was not appended at the final active slot");
        }
        int beforeDim = beforeCount * stateCount + 1;
        adjointAfter[beforeDim - 1] = adjointAfter[afterDim - 1];
    }

    /**
     * {@code afterDim} is the logical size of (i.e. number of meaningful leading
     * elements in) {@code adjointAfter} -- it cannot be read off {@code
     * adjointAfter.length}, since that array is a reused, growth-only buffer
     * from {@link #reverse}'s operation-level ping-pong and may be physically
     * larger than the current logical dimension.
     */
    private void reverseCoalescentInto(MascotCore.CoalescentTape tape, double[] adjointAfter, int afterDim,
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

            // Adjoint node-state score: pi_s = p_s * (centered + ellAdjoint), the
            // exact derivative of logL with respect to a hypothetical local
            // log-weight on this coalescent event's parent state s (see
            // MASCOT_ADJOINT_ANCESTRAL_RECONSTRUCTION.md Section 3 for the
            // derivation). It both propagates the child adjoints and
            // accumulates the population-size gradient below -- previously
            // computed as two separate terms ("bar" plus an ellAdjoint term)
            // that summed to the same value; now computed once and reused.
            double nodeStateScore = tape.parentProbabilities[s] * (centered + ellAdjoint);

            // First write to each child slot must assign (=), not accumulate (+=):
            // reverseCoalescentInto may write into a reused buffer (see reverse()'s
            // operation-level ping-pong) whose child1Index/child2Index positions the
            // compaction loop above never touches, so they can hold stale data from an
            // earlier reverse-pass operation.
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
