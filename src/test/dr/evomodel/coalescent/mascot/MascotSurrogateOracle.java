/*
 * MascotSurrogateOracle.java
 *
 * Copyright (c) 2002-2026 Alexei Drummond, Andrew Rambaut and Marc Suchard
 */

package test.dr.evomodel.coalescent.mascot;

import dr.evolution.tree.NodeRef;
import dr.evomodel.bigfasttree.BigFastTreeIntervals;
import dr.evomodel.coalescent.StructuredCoalescentSchedule;
import dr.evomodel.tree.TreeModel;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Independent reference for Pi_h, the frozen-environment MASCOT surrogate on MASCOT's RK4 grid that
 * {@code MascotStochasticMapper} samples. It shares only the tree-interval schedule (event order and interval
 * lengths) with production code. The RK4 forward pass (with per-step renormalization), the per-cell environment
 * lambda_i = invN_e sum_{j != i} (y0 + 2 y2 + 2 y3 + y4)_j / 6, and every Pi_h quantity are re-implemented here, using
 * Pade matrix exponentials and Van Loan block exponentials instead of the production uniformization series.
 * <p>
 * Conventions: a branch is identified by its lower node; jump counts use backward-time labels, entry a K + b counting
 * changes from a (below) to b (above).
 */
final class MascotSurrogateOracle {

    private static final double TIME_TOLERANCE = 1.0e-14;

    final int K;
    final int nodeCount;
    final int root;
    final double logLikelihood;

    private final double[] heights;
    private final int[] parent;
    private final int[] coalescentEpoch;
    private final List<int[]> coalescences = new ArrayList<>();   // {parent, child1, child2} in schedule order
    private final double[][] nodeFilter;
    private final double[][] nodeBeta;
    private final List<List<Cell>> cells;
    private final double[][] inversePopulation;

    private static final class Cell {
        double start;
        double end;
        double length;
        double[][] generator;
        double[][] transition;
        double[] alphaStart;
        double[] betaTop;
    }

    /**
     * @param tipPartials   by node number (tips only); normalized here
     * @param epochTimes    grid points (heights) separating epochs; empty for one epoch
     * @param migration     [epoch][K K] backward rates source -> sink; diagonal ignored
     * @param popSizes      [epoch][K] natural-scale population sizes
     * @param clockRates    branch rate above each node by node number; null for 1
     */
    MascotSurrogateOracle(TreeModel tree, double[][] tipPartials, int K, double[] epochTimes, double[][] migration,
                          double[][] popSizes, double[] clockRates, double maxStep) {
        this.K = K;
        this.nodeCount = tree.getNodeCount();
        int epochCount = epochTimes.length + 1;
        double[] boundaries = new double[epochCount + 1];
        for (int e = 1; e < epochCount; e++) {
            boundaries[e] = epochTimes[e - 1];
        }
        boundaries[epochCount] = Double.POSITIVE_INFINITY;

        double[][] M = new double[epochCount][K * K];
        inversePopulation = new double[epochCount][K];
        for (int e = 0; e < epochCount; e++) {
            for (int a = 0; a < K; a++) {
                double rowSum = 0.0;
                for (int b = 0; b < K; b++) {
                    if (a != b) {
                        M[e][a * K + b] = migration[e][a * K + b];
                        rowSum += migration[e][a * K + b];
                    }
                }
                M[e][a * K + a] = -rowSum;
                inversePopulation[e][a] = 1.0 / popSizes[e][a];
            }
        }

        heights = new double[nodeCount];
        parent = new int[nodeCount];
        Arrays.fill(parent, -1);
        coalescentEpoch = new int[nodeCount];
        for (int i = 0; i < nodeCount; i++) {
            NodeRef node = tree.getNode(i);
            heights[node.getNumber()] = tree.getNodeHeight(node);
        }
        nodeFilter = new double[nodeCount][];
        nodeBeta = new double[nodeCount][];
        cells = new ArrayList<>();
        for (int i = 0; i < nodeCount; i++) {
            cells.add(new ArrayList<>());
        }

        StructuredCoalescentSchedule schedule = StructuredCoalescentSchedule.fromTreeIntervals(tree,
                new BigFastTreeIntervals(tree), true, false);

        // MASCOT forward state (p per active lineage, in insertion order) and the Pi_h filter alpha per lineage.
        List<Integer> active = new ArrayList<>();
        double[][] p = new double[nodeCount][];
        double[][] alpha = new double[nodeCount][];
        double ell = 0.0;

        double currentTime = schedule.initialTime;
        addSample(schedule.initialSampleNode, tipPartials, active, p, alpha);
        for (int interval = 0; interval < schedule.intervalLengths.length; interval++) {
            double intervalEnd = currentTime + schedule.intervalLengths[interval];
            while (intervalEnd > currentTime + TIME_TOLERANCE) {
                double segmentEnd = nextBoundaryAfter(boundaries, currentTime, intervalEnd);
                int epoch = epochAt(boundaries, currentTime + TIME_TOLERANCE);
                ell = integrateSegment(active, p, alpha, ell, currentTime, segmentEnd, M[epoch],
                        inversePopulation[epoch], clockRates, maxStep);
                currentTime = segmentEnd;
            }
            currentTime = intervalEnd;

            if (schedule.eventTypes[interval] == StructuredCoalescentSchedule.SAMPLE) {
                addSample(schedule.sampleNodes[interval], tipPartials, active, p, alpha);
            } else {
                int c1 = schedule.coalescentChild1[interval];
                int c2 = schedule.coalescentChild2[interval];
                int v = schedule.coalescentParent[interval];
                int epoch = epochAt(boundaries, currentTime);
                double[] q = inversePopulation[epoch];
                double lambda = 0.0;
                double[] pv = new double[K];
                double[] av = new double[K];
                for (int s = 0; s < K; s++) {
                    pv[s] = p[c1][s] * p[c2][s] * q[s];
                    lambda += pv[s];
                    av[s] = alpha[c1][s] * alpha[c2][s] * q[s];
                }
                for (int s = 0; s < K; s++) {
                    pv[s] /= lambda;
                }
                ell += Math.log(lambda);
                normalize(av);
                p[v] = pv;
                alpha[v] = av;
                nodeFilter[v] = av.clone();
                parent[c1] = v;
                parent[c2] = v;
                coalescentEpoch[v] = epoch;
                coalescences.add(new int[]{v, c1, c2});
                active.remove(Integer.valueOf(c1));
                active.remove(Integer.valueOf(c2));
                active.add(v);
            }
        }
        this.logLikelihood = ell;
        this.root = coalescences.get(coalescences.size() - 1)[0];

        snapCells();
        backwardPass();
    }

    // ---------------------------------------------------------------------------------------------------------
    // Pi_h quantities
    // ---------------------------------------------------------------------------------------------------------

    /** P(node = s), for every node including tips; [node][K]. */
    double[][] nodeMarginals() {
        double[][] out = new double[nodeCount][K];
        for (int v = 0; v < nodeCount; v++) {
            for (int s = 0; s < K; s++) {
                out[v][s] = nodeFilter[v][s] * nodeBeta[v][s];
            }
            normalize(out[v]);
        }
        return out;
    }

    /** E[number of backward-time jumps a -> b on the branch above node]; [node][K K], null for the root. */
    double[][] expectedJumps() {
        double[][] out = new double[nodeCount][];
        for (int v = 0; v < nodeCount; v++) {
            if (v == root) {
                continue;
            }
            out[v] = new double[K * K];
            for (Cell cell : cells.get(v)) {
                double denominator = bilinear(cell.alphaStart, cell.transition, cell.betaTop);
                for (int a = 0; a < K; a++) {
                    for (int b = 0; b < K; b++) {
                        if (a == b || cell.generator[a][b] == 0.0) {
                            continue;
                        }
                        double[][] register = new double[K][K];
                        register[a][b] = cell.generator[a][b];
                        double[][] integral = vanLoan(cell.generator, register, cell.length);
                        out[v][a * K + b] += bilinear(cell.alphaStart, integral, cell.betaTop) / denominator;
                    }
                }
            }
        }
        return out;
    }

    /** E[time spent in each state on the branch above node]; [node][K], null for the root. */
    double[][] expectedDwell() {
        double[][] out = new double[nodeCount][];
        for (int v = 0; v < nodeCount; v++) {
            if (v == root) {
                continue;
            }
            out[v] = new double[K];
            for (Cell cell : cells.get(v)) {
                double denominator = bilinear(cell.alphaStart, cell.transition, cell.betaTop);
                for (int a = 0; a < K; a++) {
                    double[][] register = new double[K][K];
                    register[a][a] = 1.0;
                    double[][] integral = vanLoan(cell.generator, register, cell.length);
                    out[v][a] += bilinear(cell.alphaStart, integral, cell.betaTop) / denominator;
                }
            }
        }
        return out;
    }

    /** E[z^N] for the total number of jumps N on the branch above node; 1 for the root. */
    double[] jumpCountGeneratingFunction(double z) {
        double[] out = new double[nodeCount];
        for (int v = 0; v < nodeCount; v++) {
            List<Cell> branch = cells.get(v);
            if (v == root || branch.isEmpty()) {
                out[v] = 1.0;
                continue;
            }
            double[] weighted = nodeFilter[v].clone();
            double[] plain = nodeFilter[v].clone();
            for (Cell cell : branch) {
                double[][] scaled = new double[K][K];
                for (int a = 0; a < K; a++) {
                    for (int b = 0; b < K; b++) {
                        scaled[a][b] = (a == b ? 1.0 : z) * cell.generator[a][b] * cell.length;
                    }
                }
                weighted = rowTimes(weighted, MatrixExp.expm(scaled));
                plain = rowTimes(plain, cell.transition);
                double scale = sum(plain);
                for (int s = 0; s < K; s++) {
                    weighted[s] /= scale;
                    plain[s] /= scale;
                }
            }
            double[] top = branch.get(branch.size() - 1).betaTop;
            out[v] = dot(weighted, top) / dot(plain, top);
        }
        return out;
    }

    /**
     * E[number of jumps on the branch above node whose offset within their cell, (t - start) / (end - start), is
     * below {@code fraction}]; 0 for the root.
     */
    double[] expectedEarlyJumps(double fraction) {
        double[] out = new double[nodeCount];
        for (int v = 0; v < nodeCount; v++) {
            if (v == root) {
                continue;
            }
            for (Cell cell : cells.get(v)) {
                double[][] offDiagonal = new double[K][K];
                for (int a = 0; a < K; a++) {
                    for (int b = 0; b < K; b++) {
                        if (a != b) {
                            offDiagonal[a][b] = cell.generator[a][b];
                        }
                    }
                }
                double early = fraction * cell.length;
                double[][] integral = vanLoan(cell.generator, offDiagonal, early);
                double[][] rest = MatrixExp.expm(scaled(cell.generator, cell.length - early));
                double numerator = bilinear(cell.alphaStart, MatrixExp.mul(integral, rest), cell.betaTop);
                out[v] += numerator / bilinear(cell.alphaStart, cell.transition, cell.betaTop);
            }
        }
        return out;
    }

    /** Joint law of the internal nodes listed in {@code order}; index = sum_q s_q K^q. */
    double[] internalJoint(int[] order) {
        double[][][] branchMatrix = new double[nodeCount][][];
        for (int v = 0; v < nodeCount; v++) {
            if (v == root) {
                continue;
            }
            double[][] product = MatrixExp.identity(K);
            for (Cell cell : cells.get(v)) {
                product = MatrixExp.mul(product, cell.transition);
            }
            branchMatrix[v] = product;
        }
        int[] position = new int[nodeCount];
        Arrays.fill(position, -1);
        for (int q = 0; q < order.length; q++) {
            position[order[q]] = q;
        }
        int size = 1;
        for (int q = 0; q < order.length; q++) {
            size *= K;
        }
        double[] out = new double[size];
        int[] state = new int[nodeCount];
        double total = 0.0;
        for (int x = 0; x < size; x++) {
            int remainder = x;
            for (int q = 0; q < order.length; q++) {
                state[order[q]] = remainder % K;
                remainder /= K;
            }
            double weight = 1.0;
            for (int[] coalescence : coalescences) {
                int v = coalescence[0];
                weight *= inversePopulation[coalescentEpoch[v]][state[v]];
                for (int c = 1; c <= 2; c++) {
                    int child = coalescence[c];
                    if (position[child] >= 0) {
                        weight *= branchMatrix[child][state[child]][state[v]];
                    } else {
                        double tip = 0.0;
                        for (int s = 0; s < K; s++) {
                            tip += nodeFilter[child][s] * branchMatrix[child][s][state[v]];
                        }
                        weight *= tip;
                    }
                }
            }
            out[x] = weight;
            total += weight;
        }
        for (int x = 0; x < size; x++) {
            out[x] /= total;
        }
        return out;
    }

    /** Cell boundaries of the branch above node: {starts, ends}. */
    double[][] cellBoundaries(int node) {
        List<Cell> branch = cells.get(node);
        double[] starts = new double[branch.size()];
        double[] ends = new double[branch.size()];
        for (int i = 0; i < branch.size(); i++) {
            starts[i] = branch.get(i).start;
            ends[i] = branch.get(i).end;
        }
        return new double[][]{starts, ends};
    }

    int getParent(int node) {
        return parent[node];
    }

    double getHeight(int node) {
        return heights[node];
    }

    // ---------------------------------------------------------------------------------------------------------
    // Forward pass
    // ---------------------------------------------------------------------------------------------------------

    private void addSample(int node, double[][] tipPartials, List<Integer> active, double[][] p, double[][] alpha) {
        double[] partial = tipPartials[node].clone();
        normalize(partial);
        p[node] = partial.clone();
        alpha[node] = partial.clone();
        nodeFilter[node] = partial.clone();
        active.add(node);
    }

    private double integrateSegment(List<Integer> active, double[][] p, double[][] alpha, double ell, double start,
                                    double end, double[] M, double[] q, double[] clockRates, double maxStep) {
        int steps = Math.max(1, (int) Math.ceil((end - start) / maxStep));
        double h = (end - start) / steps;
        int n = active.size();
        double[] c = new double[n];
        double[] y = new double[n * K + 1];
        for (int i = 0; i < n; i++) {
            int lineage = active.get(i);
            c[i] = clockRates == null ? 1.0 : clockRates[lineage];
            System.arraycopy(p[lineage], 0, y, i * K, K);
        }
        y[n * K] = ell;

        double[] k1 = new double[y.length];
        double[] k2 = new double[y.length];
        double[] k3 = new double[y.length];
        double[] k4 = new double[y.length];
        double[] y2 = new double[y.length];
        double[] y3 = new double[y.length];
        double[] y4 = new double[y.length];
        for (int step = 0; step < steps; step++) {
            rhs(y, n, c, M, q, k1);
            for (int i = 0; i < y.length; i++) {
                y2[i] = y[i] + 0.5 * h * k1[i];
            }
            rhs(y2, n, c, M, q, k2);
            for (int i = 0; i < y.length; i++) {
                y3[i] = y[i] + 0.5 * h * k2[i];
            }
            rhs(y3, n, c, M, q, k3);
            for (int i = 0; i < y.length; i++) {
                y4[i] = y[i] + h * k3[i];
            }
            rhs(y4, n, c, M, q, k4);

            double[] sums = new double[K];
            double[][] averages = new double[n][K];
            for (int i = 0; i < n; i++) {
                for (int s = 0; s < K; s++) {
                    int index = i * K + s;
                    averages[i][s] = (y[index] + 2.0 * y2[index] + 2.0 * y3[index] + y4[index]) / 6.0;
                    sums[s] += averages[i][s];
                }
            }
            double cellStart = start + step * h;
            for (int i = 0; i < n; i++) {
                int lineage = active.get(i);
                Cell cell = new Cell();
                cell.start = cellStart;
                cell.length = h;
                cell.generator = new double[K][K];
                for (int a = 0; a < K; a++) {
                    for (int b = 0; b < K; b++) {
                        cell.generator[a][b] = c[i] * M[a * K + b];
                    }
                    cell.generator[a][a] -= Math.max(0.0, (sums[a] - averages[i][a]) * q[a]);
                }
                cell.transition = MatrixExp.expm(scaled(cell.generator, h));
                cell.alphaStart = alpha[lineage].clone();
                alpha[lineage] = rowTimes(alpha[lineage], cell.transition);
                normalize(alpha[lineage]);
                cells.get(lineage).add(cell);
            }

            for (int i = 0; i < y.length; i++) {
                y[i] = y[i] + (h / 6.0) * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
            }
            for (int i = 0; i < n; i++) {
                double total = 0.0;
                for (int s = 0; s < K; s++) {
                    total += y[i * K + s];
                }
                for (int s = 0; s < K; s++) {
                    y[i * K + s] /= total;
                }
            }
        }
        for (int i = 0; i < n; i++) {
            System.arraycopy(y, i * K, p[active.get(i)], 0, K);
        }
        return y[n * K];
    }

    /** MASCOT right-hand side: migration, conditioning on no coalescence, and the log-likelihood hazard. */
    private void rhs(double[] y, int n, double[] c, double[] M, double[] q, double[] out) {
        double[] sums = new double[K];
        double[] squares = new double[K];
        for (int i = 0; i < n; i++) {
            for (int s = 0; s < K; s++) {
                double v = y[i * K + s];
                sums[s] += v;
                squares[s] += v * v;
            }
        }
        for (int i = 0; i < n; i++) {
            double r = 0.0;
            double[] hazard = new double[K];
            for (int s = 0; s < K; s++) {
                hazard[s] = (sums[s] - y[i * K + s]) * q[s];
                r += y[i * K + s] * hazard[s];
            }
            for (int b = 0; b < K; b++) {
                double migration = 0.0;
                for (int a = 0; a < K; a++) {
                    migration += y[i * K + a] * M[a * K + b];
                }
                out[i * K + b] = c[i] * migration - y[i * K + b] * (hazard[b] - r);
            }
        }
        double total = 0.0;
        for (int s = 0; s < K; s++) {
            total += 0.5 * q[s] * (sums[s] * sums[s] - squares[s]);
        }
        out[n * K] = -total;
    }

    private static double nextBoundaryAfter(double[] boundaries, double start, double stop) {
        for (int i = 1; i < boundaries.length; i++) {
            if (boundaries[i] > start + TIME_TOLERANCE) {
                return boundaries[i] <= stop + TIME_TOLERANCE ? Math.min(boundaries[i], stop) : stop;
            }
        }
        return stop;
    }

    private static int epochAt(double[] boundaries, double t) {
        int epoch = 0;
        while (epoch + 2 < boundaries.length && t >= boundaries[epoch + 1] - TIME_TOLERANCE) {
            epoch++;
        }
        return epoch;
    }

    /** Snaps each branch's cells to the tree: lowest starts at the node, each ends where the next starts. */
    private void snapCells() {
        for (int v = 0; v < nodeCount; v++) {
            List<Cell> branch = cells.get(v);
            if (branch.isEmpty()) {
                continue;
            }
            branch.get(0).start = heights[v];
            for (int i = 0; i + 1 < branch.size(); i++) {
                branch.get(i).end = branch.get(i + 1).start;
            }
            branch.get(branch.size() - 1).end = heights[parent[v]];
        }
    }

    private void backwardPass() {
        nodeBeta[root] = ones();
        for (int i = coalescences.size() - 1; i >= 0; i--) {
            int[] coalescence = coalescences.get(i);
            int v = coalescence[0];
            for (int c = 1; c <= 2; c++) {
                int child = coalescence[c];
                int sibling = coalescence[3 - c];
                double[] siblingTop = topFilter(sibling);
                double[] beta = new double[K];
                for (int s = 0; s < K; s++) {
                    beta[s] = nodeBeta[v][s] * siblingTop[s] * inversePopulation[coalescentEpoch[v]][s];
                }
                normalize(beta);
                List<Cell> branch = cells.get(child);
                for (int k = branch.size() - 1; k >= 0; k--) {
                    Cell cell = branch.get(k);
                    cell.betaTop = beta.clone();
                    beta = columnTimes(cell.transition, beta);
                    normalize(beta);
                }
                nodeBeta[child] = beta;
            }
        }
    }

    /** Pi_h filter at the top of the branch above node. */
    private double[] topFilter(int node) {
        double[] filter = nodeFilter[node].clone();
        for (Cell cell : cells.get(node)) {
            filter = rowTimes(filter, cell.transition);
            normalize(filter);
        }
        return filter;
    }

    // ---------------------------------------------------------------------------------------------------------
    // Linear algebra
    // ---------------------------------------------------------------------------------------------------------

    /** int_0^t exp(G u) B exp(G (t - u)) du, as the top-right block of exp([[G, B], [0, G]] t). */
    private double[][] vanLoan(double[][] G, double[][] B, double t) {
        double[][] block = new double[2 * K][2 * K];
        for (int a = 0; a < K; a++) {
            for (int b = 0; b < K; b++) {
                block[a][b] = G[a][b] * t;
                block[K + a][K + b] = G[a][b] * t;
                block[a][K + b] = B[a][b] * t;
            }
        }
        double[][] exponential = MatrixExp.expm(block);
        double[][] out = new double[K][K];
        for (int a = 0; a < K; a++) {
            for (int b = 0; b < K; b++) {
                out[a][b] = exponential[a][K + b];
            }
        }
        return out;
    }

    private double bilinear(double[] row, double[][] matrix, double[] column) {
        double total = 0.0;
        for (int a = 0; a < K; a++) {
            for (int b = 0; b < K; b++) {
                total += row[a] * matrix[a][b] * column[b];
            }
        }
        return total;
    }

    private double[] rowTimes(double[] row, double[][] matrix) {
        double[] out = new double[K];
        for (int b = 0; b < K; b++) {
            for (int a = 0; a < K; a++) {
                out[b] += row[a] * matrix[a][b];
            }
        }
        return out;
    }

    private double[] columnTimes(double[][] matrix, double[] column) {
        double[] out = new double[K];
        for (int a = 0; a < K; a++) {
            for (int b = 0; b < K; b++) {
                out[a] += matrix[a][b] * column[b];
            }
        }
        return out;
    }

    private static double[][] scaled(double[][] matrix, double factor) {
        double[][] out = new double[matrix.length][matrix.length];
        for (int a = 0; a < matrix.length; a++) {
            for (int b = 0; b < matrix.length; b++) {
                out[a][b] = matrix[a][b] * factor;
            }
        }
        return out;
    }

    private double[] ones() {
        double[] out = new double[K];
        Arrays.fill(out, 1.0);
        return out;
    }

    private static double sum(double[] values) {
        double total = 0.0;
        for (double v : values) {
            total += v;
        }
        return total;
    }

    private static double dot(double[] a, double[] b) {
        double total = 0.0;
        for (int i = 0; i < a.length; i++) {
            total += a[i] * b[i];
        }
        return total;
    }

    private static void normalize(double[] values) {
        double total = sum(values);
        if (!(total > 0.0) || !Double.isFinite(total)) {
            throw new IllegalStateException("cannot normalize a vector with total " + total);
        }
        for (int i = 0; i < values.length; i++) {
            values[i] /= total;
        }
    }
}
