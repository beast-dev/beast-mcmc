package test.dr.evomodel.coalescent.mascot;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Brute-force EXACT structured coalescent on a fixed, time-scaled, binary tree with labelled lineages.
 * No BEAST dependencies.
 *
 * Conventions (chosen to match MASCOT's GenericMascotLikelihoodDelegate / MascotForwardEvaluator):
 *  - time t = node height, measured backward from the present (tips at t >= 0, root at the top).
 *  - M[e][a*K+b] (a != b) = backward-in-time migration rate of ONE lineage from deme a ("source") to deme b
 *    ("sink") in epoch e; diagonal entries are ignored (row sums are rebuilt).  This is exactly MASCOT's
 *    EpochRates.migrationMatrix and the row-major off-diagonal order of a MASCOT migration-rate Parameter.
 *  - lineage i's migration rates are multiplied by clockRate[i] (branch rate of the branch above node i).
 *  - coalescence rate of a pair of lineages that are both in deme s = 1 / N[e][s]   (N natural scale).
 *  - epochStarts = {0, g_1, g_2, ...}: epoch e covers [epochStarts[e], epochStarts[e+1]).
 *  - tip partials are normalised to sum to 1 (as StructuredCoalescentTipData does), so an ambiguous tip is a
 *    uniform mixture over its allowed demes.
 *
 * Joint process: the vector f over K^m configurations of the m active lineages solves df/dt = f G with
 *   G = sum_i c_i (M acting on coordinate i) - diag( sum_{i<j, s_i = s_j = s} 1/N_s )
 * between events; at a coalescence of i,j into p: f'(.., s_p = s) = f(.., s_i = s_j = s) / N_s;
 * at a sample event the new lineage's coordinate is multiplied in with its tip partial; Z = sum_s f_root(s).
 * logZ is the exact log density of the genealogy given the tip locations.
 *
 * exp(G dt) is applied to vectors by uniformisation (all terms non-negative, no cancellation), split into
 * sub-steps with lambda*dt <= 40, truncated when the Poisson tail is < 1e-22.
 */
public final class ExactStructuredCoalescent {

    // ------------------------------------------------------------------------------------------------
    // Inputs
    // ------------------------------------------------------------------------------------------------

    public static final class Model {
        public final int K;
        public final double[] epochStarts;  // length E, epochStarts[0] == 0
        public final double[][] M;          // [E][K*K], off-diagonal backward-time rates source->sink
        public final double[][] N;          // [E][K], natural-scale population sizes

        public Model(int K, double[] epochStarts, double[][] M, double[][] N) {
            if (epochStarts.length < 1 || epochStarts[0] != 0.0) {
                throw new IllegalArgumentException("epochStarts must start with 0");
            }
            if (M.length != epochStarts.length || N.length != epochStarts.length) {
                throw new IllegalArgumentException("need one M and one N block per epoch");
            }
            for (int e = 1; e < epochStarts.length; e++) {
                if (!(epochStarts[e] > epochStarts[e - 1])) throw new IllegalArgumentException("epochs not increasing");
            }
            this.K = K;
            this.epochStarts = epochStarts.clone();
            this.M = new double[M.length][];
            this.N = new double[N.length][];
            for (int e = 0; e < M.length; e++) {
                if (M[e].length != K * K || N[e].length != K) throw new IllegalArgumentException("bad block size");
                this.M[e] = M[e].clone();
                this.N[e] = N[e].clone();
                for (int a = 0; a < K; a++) {
                    this.M[e][a * K + a] = 0.0;
                    if (!(N[e][a] > 0)) throw new IllegalArgumentException("N must be positive");
                    for (int b = 0; b < K; b++) {
                        if (a != b && !(M[e][a * K + b] >= 0)) throw new IllegalArgumentException("negative rate");
                    }
                }
            }
        }

        /** Single epoch convenience constructor. */
        public static Model constant(int K, double[] M, double[] N) {
            return new Model(K, new double[]{0.0}, new double[][]{M}, new double[][]{N});
        }

        public int epochCount() {
            return epochStarts.length;
        }

        public int epochAt(double t) {
            int e = 0;
            while (e + 1 < epochStarts.length && t >= epochStarts[e + 1]) e++;
            return e;
        }

        /** Model with every migration rate multiplied by s and every N multiplied by n (for sweeps). */
        public Model scaled(double s, double n) {
            double[][] m2 = new double[M.length][];
            double[][] n2 = new double[N.length][];
            for (int e = 0; e < M.length; e++) {
                m2[e] = M[e].clone();
                n2[e] = N[e].clone();
                for (int i = 0; i < m2[e].length; i++) m2[e][i] *= s;
                for (int i = 0; i < n2[e].length; i++) n2[e][i] *= n;
            }
            return new Model(K, epochStarts, m2, n2);
        }
    }

    public static final class TreeSpec {
        public final int nodeCount;
        public final int[] parent;       // -1 for the root
        public final int[] child1, child2; // -1 for tips
        public final double[] height;
        public final double[][] tipPartial; // normalised; null for internal nodes
        public final double[] clockRate;  // rate of the branch above each node (ignored for the root)
        public final String[] label;      // taxon id for tips (may be null for internal nodes)
        public final int root;

        public TreeSpec(int[] parent, int[] child1, int[] child2, double[] height, double[][] tipPartial,
                        double[] clockRate, String[] label) {
            this.nodeCount = parent.length;
            this.parent = parent.clone();
            this.child1 = child1.clone();
            this.child2 = child2.clone();
            this.height = height.clone();
            this.tipPartial = new double[nodeCount][];
            this.clockRate = clockRate == null ? filled(nodeCount, 1.0) : clockRate.clone();
            this.label = label == null ? new String[nodeCount] : label.clone();
            int r = -1;
            for (int i = 0; i < nodeCount; i++) {
                if (parent[i] < 0) {
                    if (r >= 0) throw new IllegalArgumentException("two roots");
                    r = i;
                }
                boolean tip = child1[i] < 0;
                if (tip != (child2[i] < 0)) throw new IllegalArgumentException("non-binary node " + i);
                if (tip) {
                    double[] w = tipPartial[i];
                    if (w == null) throw new IllegalArgumentException("missing tip partial for node " + i);
                    double s = 0;
                    for (double v : w) {
                        if (!(v >= 0)) throw new IllegalArgumentException("bad tip partial");
                        s += v;
                    }
                    this.tipPartial[i] = new double[w.length];
                    for (int k = 0; k < w.length; k++) this.tipPartial[i][k] = w[k] / s;
                } else if (parent[child1[i]] != i || parent[child2[i]] != i) {
                    throw new IllegalArgumentException("inconsistent parent/child at " + i);
                }
            }
            this.root = r;
        }

        public boolean isTip(int i) {
            return child1[i] < 0;
        }

        public List<Integer> internalNodes() {
            List<Integer> out = new ArrayList<>();
            for (int i = 0; i < nodeCount; i++) if (!isTip(i)) out.add(i);
            return out;
        }

        /** Sorted, comma-separated tip labels below node i (the clade key). */
        public String cladeKey(int i) {
            List<String> names = new ArrayList<>();
            collect(i, names);
            java.util.Collections.sort(names);
            return String.join(",", names);
        }

        private void collect(int i, List<String> names) {
            if (isTip(i)) {
                names.add(label[i] == null ? ("#" + i) : label[i]);
            } else {
                collect(child1[i], names);
                collect(child2[i], names);
            }
        }
    }

    /**
     * Optional forward-mode perturbation: the dual channel carries dZ/dtheta where
     *  - JUMP: theta multiplies only the off-diagonal rate source->sink of lineage `node` (all epochs);
     *          E[#jumps source->sink on branch above node | data] = (dZ/dtheta)/Z at theta=1.
     *  - DWELL: theta is an extra killing rate on lineage `node` while it is in deme `source`;
     *          E[time in source on branch above node | data] = -(dZ/dtheta)/Z at theta=0.
     */
    public static final class Perturbation {
        public static final int JUMP = 0, DWELL = 1;
        final int type, node, source, sink;

        public Perturbation(int type, int node, int source, int sink) {
            this.type = type;
            this.node = node;
            this.source = source;
            this.sink = sink;
        }
    }

    // ------------------------------------------------------------------------------------------------
    // Engine
    // ------------------------------------------------------------------------------------------------

    private final Model model;
    private final TreeSpec tree;
    private final int K;
    private final int[] order; // nodes sorted by height (tips before internals on ties)

    /** Tolerance for treating two event heights as simultaneous. */
    public static final double TIME_TOL = 1e-12;
    static final double MAX_MU = 40.0;
    static final double TAIL = 1e-22;

    public ExactStructuredCoalescent(Model model, TreeSpec tree) {
        this.model = model;
        this.tree = tree;
        this.K = model.K;
        Integer[] idx = new Integer[tree.nodeCount];
        for (int i = 0; i < idx.length; i++) idx[i] = i;
        Arrays.sort(idx, (a, b) -> {
            int c = Double.compare(tree.height[a], tree.height[b]);
            if (c != 0 && Math.abs(tree.height[a] - tree.height[b]) > TIME_TOL) return c;
            boolean ta = tree.isTip(a), tb = tree.isTip(b);
            if (ta != tb) return ta ? -1 : 1;
            return c;
        });
        order = new int[idx.length];
        for (int i = 0; i < idx.length; i++) order[i] = idx[i];
        // sanity: no coalescence simultaneous with another event
        for (int i = 1; i < order.length; i++) {
            int a = order[i - 1], b = order[i];
            if (Math.abs(tree.height[a] - tree.height[b]) <= TIME_TOL && !(tree.isTip(a) && tree.isTip(b))) {
                throw new IllegalArgumentException("simultaneous events involving a coalescence at nodes " + a + "," + b);
            }
        }
        for (int i = 0; i < tree.nodeCount; i++) {
            if (tree.isTip(i) && tree.tipPartial[i].length != K) throw new IllegalArgumentException("tip partial size");
            if (!tree.isTip(i)) {
                if (!(tree.height[i] > tree.height[tree.child1[i]] && tree.height[i] > tree.height[tree.child2[i]])) {
                    throw new IllegalArgumentException("parent not older than child at node " + i);
                }
            }
        }
    }

    /** Result of one forward pass. */
    public static final class Pass {
        public final double logZ;      // log of Z (primal)
        public final double dualRatio; // (dZ/dtheta) / Z  (0 if no perturbation)

        Pass(double logZ, double dualRatio) {
            this.logZ = logZ;
            this.dualRatio = dualRatio;
        }
    }

    public double logLikelihood() {
        return forward(null, null).logZ;
    }

    /**
     * One forward pass (tips -> root).
     * @param nodeWeights optional [node][K] non-negative multipliers applied to the parent density at each
     *                    coalescence (an indicator e_s gives Z restricted to node = s).
     * @param pert        optional dual-number perturbation (see {@link Perturbation}).
     */
    public Pass forward(double[][] nodeWeights, Perturbation pert) {
        int[] active = new int[0];
        double[] f = new double[]{1.0};
        double[] g = pert == null ? null : new double[]{0.0};
        double logScale = 0.0;
        double t = tree.height[order[0]];

        for (int oi = 0; oi < order.length; oi++) {
            int node = order[oi];
            double tNext = tree.height[node];
            if (tNext > t + TIME_TOL) {
                // propagate across [t, tNext], split at epoch boundaries
                double a = t;
                while (a < tNext - TIME_TOL) {
                    int e = model.epochAt(a + TIME_TOL);
                    double b = tNext;
                    if (e + 1 < model.epochCount() && model.epochStarts[e + 1] < b) b = model.epochStarts[e + 1];
                    double[][] fg = propagate(active, f, g, e, b - a, pert);
                    f = fg[0];
                    g = fg[1];
                    a = b;
                }
                t = tNext;
            }
            if (tree.isTip(node)) {
                double[] w = tree.tipPartial[node];
                f = appendLineage(f, active.length, w);
                if (g != null) g = appendLineage(g, active.length, w);
                active = Arrays.copyOf(active, active.length + 1);
                active[active.length - 1] = node;
            } else {
                int e = model.epochAt(tNext);
                int i = indexOf(active, tree.child1[node]);
                int j = indexOf(active, tree.child2[node]);
                if (i < 0 || j < 0) throw new IllegalStateException("children not active at node " + node);
                double[] wt = new double[K];
                for (int s = 0; s < K; s++) {
                    wt[s] = (1.0 / model.N[e][s]) * (nodeWeights == null || nodeWeights[node] == null ? 1.0 : nodeWeights[node][s]);
                }
                f = coalesce(f, active.length, i, j, wt);
                if (g != null) g = coalesce(g, active.length, i, j, wt);
                int[] na = new int[active.length - 1];
                int k = 0;
                for (int q = 0; q < active.length; q++) if (q != i && q != j) na[k++] = active[q];
                na[k] = node;
                active = na;
            }
            // rescale
            double s = 0;
            for (double v : f) s += v;
            if (s > 0 && Double.isFinite(s)) {
                for (int q = 0; q < f.length; q++) f[q] /= s;
                if (g != null) for (int q = 0; q < g.length; q++) g[q] /= s;
                logScale += Math.log(s);
            } else if (s == 0) {
                return new Pass(Double.NEGATIVE_INFINITY, Double.NaN);
            } else {
                throw new IllegalStateException("non-finite mass " + s);
            }
        }
        if (active.length != 1) throw new IllegalStateException("ended with " + active.length + " lineages");
        double z = 0, zg = 0;
        for (int s = 0; s < K; s++) {
            z += f[s];
            if (g != null) zg += g[s];
        }
        return new Pass(logScale + Math.log(z), g == null ? 0.0 : zg / z);
    }

    // -- events --------------------------------------------------------------------------------------

    private double[] appendLineage(double[] f, int m, double[] w) {
        int size = f.length;
        double[] out = new double[size * K];
        for (int s = 0; s < K; s++) {
            int off = s * size;
            for (int x = 0; x < size; x++) out[off + x] = f[x] * w[s];
        }
        return out;
    }

    /** Remove coordinates i and j (requiring equal states) and append the parent coordinate. */
    private double[] coalesce(double[] f, int m, int i, int j, double[] wt) {
        int newM = m - 1;
        int[] pow = powers(m);
        int outSize = ipow(K, newM);
        double[] out = new double[outSize];
        int[] keep = new int[m - 2];
        int k = 0;
        for (int q = 0; q < m; q++) if (q != i && q != j) keep[k++] = q;
        int restSize = ipow(K, m - 2);
        for (int y = 0; y < restSize; y++) {
            int base = 0;
            int r = y;
            for (int q = 0; q < m - 2; q++) {
                int sq = r % K;
                r /= K;
                base += sq * pow[keep[q]];
            }
            for (int s = 0; s < K; s++) {
                int x = base + s * pow[i] + s * pow[j];
                out[y + s * restSize] = f[x] * wt[s];
            }
        }
        return out;
    }

    // -- propagation by uniformisation ------------------------------------------------------------------

    private double[][] propagate(int[] active, double[] f, double[] g, int epoch, double dt, Perturbation pert) {
        int m = active.length;
        int size = f.length;
        double[] Me = model.M[epoch];
        double[] Ne = model.N[epoch];
        double[] rowSum = new double[K];
        for (int a = 0; a < K; a++) for (int b = 0; b < K; b++) if (a != b) rowSum[a] += Me[a * K + b];
        double[] c = new double[m];
        for (int i = 0; i < m; i++) c[i] = tree.clockRate[active[i]];
        int[] pow = powers(m);

        // exit rates
        double[] r = new double[size];
        double lambda = 0;
        for (int x = 0; x < size; x++) {
            double out = 0;
            int[] s = decode(x, m);
            for (int i = 0; i < m; i++) out += c[i] * rowSum[s[i]];
            for (int i = 0; i < m; i++) for (int j = i + 1; j < m; j++) if (s[i] == s[j]) out += 1.0 / Ne[s[i]];
            r[x] = out;
            lambda = Math.max(lambda, out);
        }
        int pertLineage = -1;
        if (pert != null) pertLineage = indexOf(active, pert.node);
        if (pertLineage >= 0 && pert.type == Perturbation.DWELL) lambda += 1.0; // room for the extra diagonal
        if (lambda == 0) return new double[][]{f, g};

        int sub = Math.max(1, (int) Math.ceil(lambda * dt / MAX_MU));
        double h = dt / sub;
        double mu = lambda * h;
        for (int st = 0; st < sub; st++) {
            double[] v = f.clone();
            double[] u = g == null ? null : g.clone();
            double w = Math.exp(-mu);
            double[] accF = new double[size];
            double[] accG = g == null ? null : new double[size];
            for (int q = 0; q < size; q++) accF[q] = w * v[q];
            if (u != null) for (int q = 0; q < size; q++) accG[q] = w * u[q];
            for (int k = 1; ; k++) {
                double[] v2 = new double[size];
                double[] u2 = u == null ? null : new double[size];
                applyP(v, v2, m, c, Me, r, lambda, pow);
                if (u != null) {
                    applyP(u, u2, m, c, Me, r, lambda, pow);
                    if (pertLineage >= 0) applyPerturbationDerivative(v, u2, m, pertLineage, c, Me, lambda, pow, pert);
                }
                v = v2;
                u = u2;
                w *= mu / k;
                for (int q = 0; q < size; q++) accF[q] += w * v[q];
                if (u != null) for (int q = 0; q < size; q++) accG[q] += w * u[q];
                if (k > mu && w < TAIL * (1 + k)) break;
                if (k > 100000) throw new IllegalStateException("uniformisation did not converge");
            }
            f = accF;
            g = accG;
        }
        return new double[][]{f, g};
    }

    /** out = v P, with P = I + G/lambda. */
    private void applyP(double[] v, double[] out, int m, double[] c, double[] Me, double[] r, double lambda, int[] pow) {
        int size = v.length;
        for (int x = 0; x < size; x++) {
            double vx = v[x];
            if (vx == 0) continue;
            out[x] += vx * (1.0 - r[x] / lambda);
            int rem = x;
            for (int i = 0; i < m; i++) {
                int a = rem % K;
                rem /= K;
                if (c[i] == 0) continue;
                for (int b = 0; b < K; b++) {
                    if (b == a) continue;
                    double rate = c[i] * Me[a * K + b];
                    if (rate != 0) out[x + (b - a) * pow[i]] += vx * rate / lambda;
                }
            }
        }
    }

    /** out += v dP/dtheta for the perturbation on lineage position `li`. */
    private void applyPerturbationDerivative(double[] v, double[] out, int m, int li, double[] c, double[] Me,
                                             double lambda, int[] pow, Perturbation pert) {
        int size = v.length;
        for (int x = 0; x < size; x++) {
            int a = (x / pow[li]) % K;
            if (a != pert.source) continue;
            if (pert.type == Perturbation.JUMP) {
                double rate = c[li] * Me[pert.source * K + pert.sink];
                out[x + (pert.sink - pert.source) * pow[li]] += v[x] * rate / lambda;
            } else {
                out[x] += -v[x] / lambda;
            }
        }
    }

    // -- derived quantities --------------------------------------------------------------------------

    /** Exact posterior marginals P(node = s | tree, tips), [node][K]; tips get their (normalised) partials. */
    public double[][] marginals() {
        double logZ = logLikelihood();
        double[][] out = new double[tree.nodeCount][];
        for (int node = 0; node < tree.nodeCount; node++) {
            if (tree.isTip(node)) {
                // posterior at a tip with ambiguous data: also computable via a weight on the tip partial
                out[node] = tipPosterior(node, logZ);
                continue;
            }
            out[node] = new double[K];
            for (int s = 0; s < K; s++) {
                double[][] w = new double[tree.nodeCount][];
                w[node] = new double[K];
                w[node][s] = 1.0;
                Pass p = forward(w, null);
                out[node][s] = Math.exp(p.logZ - logZ);
            }
        }
        return out;
    }

    private double[] tipPosterior(int tip, double logZ) {
        double[] w = tree.tipPartial[tip];
        double[] out = new double[K];
        int nonzero = 0;
        for (double v : w) if (v > 0) nonzero++;
        if (nonzero <= 1) return w.clone();
        for (int s = 0; s < K; s++) {
            if (w[s] == 0) continue;
            double[][] partials = new double[tree.nodeCount][];
            for (int i = 0; i < tree.nodeCount; i++) partials[i] = tree.tipPartial[i];
            double[] e = new double[K];
            e[s] = 1.0;
            partials[tip] = e;
            TreeSpec t2 = new TreeSpec(tree.parent, tree.child1, tree.child2, tree.height, partials, tree.clockRate, tree.label);
            double lz = new ExactStructuredCoalescent(model, t2).logLikelihood();
            out[s] = w[s] * Math.exp(lz - logZ);
        }
        return out;
    }

    /**
     * Exact joint posterior over all internal-node states. Internal nodes are ordered as in
     * {@code internalOrder}; entry index = sum_q s_q K^q (first node is the fastest-varying digit).
     */
    public double[] joint(int[] internalOrder) {
        double logZ = logLikelihood();
        int n = internalOrder.length;
        int size = ipow(K, n);
        double[] out = new double[size];
        for (int x = 0; x < size; x++) {
            double[][] w = new double[tree.nodeCount][];
            int rem = x;
            for (int q = 0; q < n; q++) {
                int s = rem % K;
                rem /= K;
                w[internalOrder[q]] = new double[K];
                w[internalOrder[q]][s] = 1.0;
            }
            Pass p = forward(w, null);
            out[x] = p.logZ == Double.NEGATIVE_INFINITY ? 0.0 : Math.exp(p.logZ - logZ);
        }
        return out;
    }

    /** E[#jumps source->sink on branch above node | data]; [node][K*K] (diagonal 0, root row null). */
    public double[][] expectedJumps() {
        double[][] out = new double[tree.nodeCount][];
        for (int node = 0; node < tree.nodeCount; node++) {
            if (node == tree.root) continue;
            out[node] = new double[K * K];
            for (int a = 0; a < K; a++) for (int b = 0; b < K; b++) {
                if (a == b) continue;
                out[node][a * K + b] = forward(null, new Perturbation(Perturbation.JUMP, node, a, b)).dualRatio;
            }
        }
        return out;
    }

    /** E[time in deme a on branch above node | data]; [node][K] (root row null). Rows sum to branch length. */
    public double[][] expectedDwell() {
        double[][] out = new double[tree.nodeCount][];
        for (int node = 0; node < tree.nodeCount; node++) {
            if (node == tree.root) continue;
            out[node] = new double[K];
            for (int a = 0; a < K; a++) {
                out[node][a] = -forward(null, new Perturbation(Perturbation.DWELL, node, a, a)).dualRatio;
            }
        }
        return out;
    }

    // -- helpers -------------------------------------------------------------------------------------

    private int[] decode(int x, int m) {
        int[] s = new int[m];
        for (int i = 0; i < m; i++) {
            s[i] = x % K;
            x /= K;
        }
        return s;
    }

    private int[] powers(int m) {
        int[] p = new int[Math.max(m, 1)];
        int v = 1;
        for (int i = 0; i < m; i++) {
            p[i] = v;
            v *= K;
        }
        return p;
    }

    static int ipow(int b, int e) {
        int v = 1;
        for (int i = 0; i < e; i++) v *= b;
        return v;
    }

    private static int indexOf(int[] a, int v) {
        for (int i = 0; i < a.length; i++) if (a[i] == v) return i;
        return -1;
    }

    static double[] filled(int n, double v) {
        double[] a = new double[n];
        Arrays.fill(a, v);
        return a;
    }

    public Model getModel() {
        return model;
    }

    public TreeSpec getTree() {
        return tree;
    }

    // ------------------------------------------------------------------------------------------------
    // Dense generator for validation (single interval, given active lineages and epoch).
    // ------------------------------------------------------------------------------------------------

    /** Dense joint generator G (row-vector convention: f' = f G) for m lineages with clock rates c. */
    public static double[][] denseGenerator(int K, int m, double[] c, double[] Me, double[] Ne) {
        int size = ipow(K, m);
        double[][] G = new double[size][size];
        int[] pow = new int[m];
        int v = 1;
        for (int i = 0; i < m; i++) {
            pow[i] = v;
            v *= K;
        }
        for (int x = 0; x < size; x++) {
            int[] s = new int[m];
            int rem = x;
            for (int i = 0; i < m; i++) {
                s[i] = rem % K;
                rem /= K;
            }
            for (int i = 0; i < m; i++) {
                for (int b = 0; b < K; b++) {
                    if (b == s[i]) continue;
                    double rate = c[i] * Me[s[i] * K + b];
                    G[x][x + (b - s[i]) * pow[i]] += rate;
                    G[x][x] -= rate;
                }
            }
            for (int i = 0; i < m; i++) for (int j = i + 1; j < m; j++) if (s[i] == s[j]) G[x][x] -= 1.0 / Ne[s[i]];
        }
        return G;
    }

    /** Uniformisation action f exp(G dt) for a dense G (validation helper). */
    public static double[] uniformizeDense(double[] f, double[][] G, double dt) {
        int n = f.length;
        double lambda = 0;
        for (int i = 0; i < n; i++) lambda = Math.max(lambda, -G[i][i]);
        if (lambda == 0) return f.clone();
        int sub = Math.max(1, (int) Math.ceil(lambda * dt / MAX_MU));
        double h = dt / sub, mu = lambda * h;
        double[] cur = f.clone();
        for (int st = 0; st < sub; st++) {
            double[] v = cur.clone();
            double w = Math.exp(-mu);
            double[] acc = new double[n];
            for (int q = 0; q < n; q++) acc[q] = w * v[q];
            for (int k = 1; ; k++) {
                double[] v2 = new double[n];
                for (int x = 0; x < n; x++) {
                    if (v[x] == 0) continue;
                    for (int y = 0; y < n; y++) {
                        double p = (x == y ? 1.0 : 0.0) + G[x][y] / lambda;
                        v2[y] += v[x] * p;
                    }
                }
                v = v2;
                w *= mu / k;
                for (int q = 0; q < n; q++) acc[q] += w * v[q];
                if (k > mu && w < TAIL * (1 + k)) break;
            }
            cur = acc;
        }
        return cur;
    }
}
