/*
 * MascotHTransformMapper.java
 *
 * Copyright (c) 2002-2026 Alexei Drummond, Andrew Rambaut and Marc Suchard
 */

package dr.evomodel.coalescent.mascot;

import dr.math.MathUtils;

/**
 * Draws from the continuous frozen-environment surrogate by the Doob h-transform that BEAST2's MappedMascot uses
 * (Stolz et al. 2022, Mol Biol Evol 39:msab342). Going down a branch from its parent, a lineage in state cur jumps to
 * new at rate c M[new][cur] p(t)[new] / p(t)[cur], where p is the lineage's MASCOT probability vector: the
 * normalized filter of the surrogate, so the killing cancels and no matrix exponential is needed. The root is drawn
 * from MASCOT's root probabilities and both children of a coalescence start in the parent's state.
 * <p>
 * Between RK4 step ends p is a cubic Hermite interpolant of the taped values and of MASCOT's right-hand side there,
 * and jumps are placed by exact Lewis-Shedler thinning against a per-step bound built from the cubics' extremes.
 * Where the current state's probability gets close to zero within a step -- always at a tip whose data exclude it --
 * or thinning would need too many proposals, the rest of the step is drawn exactly from the step's Pi_h cell
 * ({@link MascotGridMapper}), started from the step's lower MASCOT vector.
 * <p>
 * The draw is exact for the surrogate defined by the interpolated p; its accuracy is limited by the interpolation and,
 * in stiff steps, by RK4 itself. It is an alternative to {@link MascotGridMapper} for comparing speed.
 */
final class MascotHTransformMapper extends MascotStochasticMapper {

    /**
     * Above this largest h mu (total exit rate times step length) the RK4 solution this sampler interpolates is too
     * coarse: at h mu = 2.84 its node frequencies were off by 0.044 while the grid sampler stayed exact for Pi_h.
     * States this stiff are drawn with {@link MascotGridMapper} instead.
     */
    static final double STIFFNESS_LIMIT = 1.0;

    /** Below this probability for the current state anywhere in the step, the step is bridged exactly instead. */
    private static final double MIN_PROBABILITY = 1.0e-6;
    /** More expected thinning proposals than this in the rest of a step also switches to the exact bridge. */
    private static final double MAX_EXPECTED_PROPOSALS = 20.0;

    // Cell payload, each block K wide: killing, Hermite coefficients c0..c3 of p(tau) = c0 + c1 tau + c2 tau^2 +
    // c3 tau^3 for tau = (t - start) / (end - start), and the minimum and maximum of each p(tau) over [0, 1].
    private static final int KILLING = 0;
    private static final int C0 = 1;
    private static final int C1 = 2;
    private static final int C2 = 3;
    private static final int C3 = 4;
    private static final int MINIMUM = 5;
    private static final int MAXIMUM = 6;
    private static final int BLOCKS = 7;

    private final double[] payload;
    private final double[] sums;
    private final double[] rates;
    private final double[] hazard;
    private double[] averages = new double[0];
    private double[] rowDerivatives = new double[0];
    private double[] rootProbabilities;

    private double[] pieceFilters = new double[0];

    MascotHTransformMapper(int stateCount) {
        super(stateCount, BLOCKS * stateCount);
        this.payload = new double[BLOCKS * stateCount];
        this.sums = new double[stateCount];
        this.rates = new double[stateCount];
        this.hazard = new double[stateCount];
        this.rootProbabilities = new double[stateCount];
    }

    @Override
    protected StochasticMapSampler sampler() {
        return StochasticMapSampler.H_TRANSFORM;
    }

    @Override
    protected void onSample(int lineage, MascotPreparedInput input) {
        // The tip partial is the first row of the lineage's first interval.
    }

    @Override
    protected void onInterval(MascotAdjointTape.Interval interval, int[] slotLineage) {
        int K = stateCount;
        int activeCount = interval.activeCount;
        int dim = interval.stateDimension;
        double h = interval.h;

        // MASCOT's right-hand side at every step end of this interval (row `steps` is the interval's final state).
        rowDerivatives = ensure(rowDerivatives, (interval.steps + 1) * activeCount * K);
        for (int row = 0; row <= interval.steps; row++) {
            rightHandSide(interval, row * dim, row * activeCount * K);
        }

        averages = ensure(averages, activeCount * K);
        for (int step = 0; step < interval.steps; step++) {
            stageAverages(interval, step, averages, sums);
            double stepStart = interval.startTime + step * h;
            for (int slot = 0; slot < activeCount; slot++) {
                double clockRate = interval.clockRates == null ? 1.0 : interval.clockRates[slot];
                stepKilling(interval, slot, sums, averages, payload, KILLING * K);
                maxRateTimesLength = Math.max(maxRateTimesLength,
                        maxExitRate(interval.epoch, clockRate, payload, KILLING * K) * h);
                int lower = step * dim + slot * K;
                int upper = (step + 1) * dim + slot * K;
                int lowerDerivative = (step * activeCount + slot) * K;
                int upperDerivative = ((step + 1) * activeCount + slot) * K;
                for (int s = 0; s < K; s++) {
                    double pA = interval.y0[lower + s];
                    double pB = interval.y0[upper + s];
                    double dA = h * rowDerivatives[lowerDerivative + s];
                    double dB = h * rowDerivatives[upperDerivative + s];
                    payload[C0 * K + s] = pA;
                    payload[C1 * K + s] = dA;
                    payload[C2 * K + s] = -3.0 * pA - 2.0 * dA + 3.0 * pB - dB;
                    payload[C3 * K + s] = 2.0 * pA + dA - 2.0 * pB + dB;
                    cubicRange(payload[C0 * K + s], payload[C1 * K + s], payload[C2 * K + s], payload[C3 * K + s],
                            payload, MINIMUM * K + s, MAXIMUM * K + s);
                }
                cells.add(slotLineage[slot], interval.epoch, stepStart, h, clockRate, payload, 0);
            }
        }
    }

    @Override
    protected void onCoalescent(int parent, int child1, int child2, MascotAdjointTape.Coalescent coalescent) {
        // The last coalescence is the root's; its normalized parent vector is MASCOT's root probabilities.
        System.arraycopy(coalescent.parentProbabilities, 0, rootProbabilities, 0, stateCount);
    }

    @Override
    protected int drawRootState() {
        return drawCategorical(rootProbabilities, 0);
    }

    @Override
    protected int drawBranch(int node, int topState, KilledUniformizedBridge.Path changes) {
        int first = cells.first(node);
        int state = topState;
        for (int cell = first + cells.countOf(node) - 1; cell >= first; cell--) {
            state = drawCell(cell, state, changes);
        }
        return state;
    }

    /** Thinning down one cell from its upper end; returns the state at its lower end. */
    private int drawCell(int cell, int state, KilledUniformizedBridge.Path changes) {
        int K = stateCount;
        double[] p = cells.payload();
        int offset = cells.payloadOffset(cell);
        int mOffset = cells.epoch(cell) * K * K;
        double clockRate = cells.clock(cell);
        double start = cells.start(cell);
        double span = cells.end(cell) - start;
        double t = cells.end(cell);

        while (true) {
            double currentMinimum = p[offset + MINIMUM * K + state];
            if (currentMinimum < MIN_PROBABILITY) {
                return bridgeRest(cell, start, t, state, changes);
            }
            double bound = 0.0;
            for (int next = 0; next < K; next++) {
                if (next != state) {
                    bound += migration[mOffset + next * K + state] * Math.max(0.0, p[offset + MAXIMUM * K + next]);
                }
            }
            bound *= clockRate / currentMinimum;
            if (bound == 0.0) {
                return state;
            }
            if (bound * (t - start) > MAX_EXPECTED_PROPOSALS) {
                return bridgeRest(cell, start, t, state, changes);
            }

            t -= MathUtils.nextExponential(bound);
            if (t <= start) {
                return state;
            }
            double tau = (t - start) / span;
            double current = cubic(p, offset, state, tau);
            double total = 0.0;
            for (int next = 0; next < K; next++) {
                rates[next] = next == state ? 0.0 :
                        clockRate * migration[mOffset + next * K + state] * Math.max(0.0, cubic(p, offset, next, tau)) /
                                current;
                total += rates[next];
            }
            if (MathUtils.nextDouble() * bound < total) {
                int next = drawCategorical(rates, 0);
                changes.add(t, state);
                state = next;
            }
        }
    }

    /**
     * Draws the state at the cell's lower end and the path up to {@code top}, where the lineage is in
     * {@code topState}, from the cell's Pi_h generator started from the cell's lower MASCOT vector; splits the
     * interval when mu h is large.
     */
    private int bridgeRest(int cell, double start, double top, int topState, KilledUniformizedBridge.Path changes) {
        int K = stateCount;
        double[] p = cells.payload();
        int offset = cells.payloadOffset(cell);
        double length = top - start;
        if (!(length > 0.0)) {
            return topState;
        }
        buildGenerator(cells.epoch(cell), cells.clock(cell), p, offset + KILLING * K);
        double mu = maxExitRate(cells.epoch(cell), cells.clock(cell), p, offset + KILLING * K);
        int pieces = Math.max(1, (int) Math.ceil(mu * length / KilledUniformizedBridge.MAX_RECOMMENDED_RATE_TIMES_LENGTH));
        double pieceLength = length / pieces;
        bridge.setGenerator(generator, pieceLength);

        pieceFilters = ensure(pieceFilters, pieces * K);
        for (int s = 0; s < K; s++) {
            pieceFilters[s] = Math.max(0.0, p[offset + C0 * K + s]);
        }
        normalize(pieceFilters, 0, "lower end of a cell");
        for (int piece = 1; piece < pieces; piece++) {
            bridge.transitionRow(pieceFilters, (piece - 1) * K, pieceFilters, piece * K);
            normalize(pieceFilters, piece * K, "inside a cell");
        }
        int state = topState;
        for (int piece = pieces - 1; piece >= 0; piece--) {
            double pieceStart = start + piece * pieceLength;
            double pieceEnd = piece == pieces - 1 ? top : start + (piece + 1) * pieceLength;
            state = bridge.sampleStartAndPath(pieceFilters, piece * K, state, pieceStart, pieceEnd, changes);
        }
        return state;
    }

    /** MASCOT's lineage right-hand side for every active slot of {@code interval} at state row {@code rowOffset}. */
    private void rightHandSide(MascotAdjointTape.Interval interval, int rowOffset, int outOffset) {
        int K = stateCount;
        int activeCount = interval.activeCount;
        int mOffset = interval.epoch * K * K;
        int qOffset = interval.epoch * K;
        double[] y = interval.y0;
        java.util.Arrays.fill(sums, 0.0);
        for (int i = 0; i < activeCount * K; i++) {
            sums[i % K] += y[rowOffset + i];
        }
        for (int slot = 0; slot < activeCount; slot++) {
            int base = rowOffset + slot * K;
            double clockRate = interval.clockRates == null ? 1.0 : interval.clockRates[slot];
            double r = 0.0;
            for (int s = 0; s < K; s++) {
                hazard[s] = (sums[s] - y[base + s]) * inversePopulation[qOffset + s];
                r += y[base + s] * hazard[s];
            }
            for (int sink = 0; sink < K; sink++) {
                double flow = 0.0;
                for (int source = 0; source < K; source++) {
                    flow += y[base + source] * migration[mOffset + source * K + sink];
                }
                rowDerivatives[outOffset + slot * K + sink] = clockRate * flow - y[base + sink] * (hazard[sink] - r);
            }
        }
    }

    private double cubic(double[] p, int offset, int state, double tau) {
        int K = stateCount;
        return p[offset + C0 * K + state] + tau * (p[offset + C1 * K + state] +
                tau * (p[offset + C2 * K + state] + tau * p[offset + C3 * K + state]));
    }

    /** Minimum and maximum of c0 + c1 x + c2 x^2 + c3 x^3 over [0, 1]. */
    private static void cubicRange(double c0, double c1, double c2, double c3, double[] out, int minIndex,
                                   int maxIndex) {
        double atZero = c0;
        double atOne = c0 + c1 + c2 + c3;
        double min = Math.min(atZero, atOne);
        double max = Math.max(atZero, atOne);
        // Critical points: c1 + 2 c2 x + 3 c3 x^2 = 0.
        double a = 3.0 * c3;
        double b = 2.0 * c2;
        double c = c1;
        double[] roots = new double[2];
        int count = 0;
        if (Math.abs(a) < 1.0e-14 * (Math.abs(b) + Math.abs(c) + 1e-300)) {
            if (b != 0.0) {
                roots[count++] = -c / b;
            }
        } else {
            double discriminant = b * b - 4.0 * a * c;
            if (discriminant >= 0.0) {
                double root = Math.sqrt(discriminant);
                double q = -0.5 * (b + Math.copySign(root, b));
                roots[count++] = q / a;
                if (q != 0.0) {
                    roots[count++] = c / q;
                }
            }
        }
        for (int i = 0; i < count; i++) {
            double x = roots[i];
            if (x > 0.0 && x < 1.0) {
                double v = c0 + x * (c1 + x * (c2 + x * c3));
                min = Math.min(min, v);
                max = Math.max(max, v);
            }
        }
        out[minIndex] = min;
        out[maxIndex] = max;
    }
}
