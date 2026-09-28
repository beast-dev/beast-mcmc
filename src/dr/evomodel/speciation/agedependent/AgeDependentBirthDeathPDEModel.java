package dr.evomodel.speciation.agedependent;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeUtils;
import dr.evomodel.speciation.agedependent.agehazard.AgeHazard;
import dr.evomodel.tree.TreeChangedEvent;
import dr.inference.model.*;
import dr.math.RungeKutta;
import dr.xml.Reportable;
import dr.util.TaskPool;
import java.util.Arrays;
import java.util.concurrent.atomic.AtomicInteger;

/**
 * @author Frederik M. Andersen
 *
 * Computes the tree likelihood for a skyline, time- and age-dependent birth-death process under
 * two speciation modes:
 *   - symmetric:  both daughter lineages reset to age 0 at speciation
 *   - asymmetric: one daughter inherits the parent's age a, the other resets to age 0
 *
 * Let L(t, a) denote the partial subtree likelihood at time t for a lineage of age a, and
 * p0(t, a) the probability that a lineage of age a at time t leaves no sampled descendants.
 * lambda(t, a) and mu(t, a) are the (skyline, age-dependent) birth and death hazards.
 *
 * Sampling is serial: tips may have height > 0 (fossils / serial samples), drawn at a time-only
 * sampling rate psi(t) (skyline), and lineages extant at the present are sampled with
 * probability rho. Setting psi = 0 and rho = 1 recovers the ultrametric process.
 *
 *   Symmetric:
 *       dL/dt  = dL/da + 2 lambda p0(t,0) L(t,0) - (lambda + mu + psi) L(t,a)
 *       dp0/dt = dp0/da + mu + lambda p0(t,0)^2 - (lambda + mu + psi) p0(t,a)
 *
 *   Asymmetric:
 *       dL/dt  = dL/da + lambda (L(t,a) p0(t,0) + p0(t,a) L(t,0)) - (lambda + mu + psi) L(t,a)
 *       dp0/dt = dp0/da + mu + lambda p0(t,a) p0(t,0) - (lambda + mu + psi) p0(t,a)
 *
 * Boundary conditions:
 *       p0(0, a) = 1 - rho
 *       L_tip(t_i, a) = psi(t_i)   for a serial sample (t_i > 0)
 *       L_tip(0, a)   = rho        for an extant sample
 *
 * Implementation note: each tip is initialized with L = 1 and the psi(t_i) / rho factors are
 * accumulated separately as logSamplingFactor in calculateLogLikelihood(). That is equivalent
 * to the boundary above — sampling factors do not interact with the PDE step, so pulling them
 * out of the per-branch L state is just bookkeeping.
 *
 * Rate form:
 *      rate(t, a) = scale(t) * h(a)
 * where h is supplied as an {@link AgeHazard} per epoch (one per epoch for both birth and death).
 *
 * Discretization: {@code maxTime} is the outer time/age bound the PDE grid is sized to
 * ({@code Na = ceil(maxTime/deltaA)}, {@code Nt = max(Na, ceil(maxTime/deltaT))}); actual grid
 * spacing may end up marginally finer than the requested deltas. {@code maxTime} additionally
 * plays the role of the process's origin time whenever {@code conditionAt == ORIGIN}, since
 * the root-branch survival term and the {@code p0} integration endpoint need an origin — when
 * {@code conditionAt == MRCA} it is used only to size the grid.
 *
 * Solver: {@code solver} selects the numerical scheme.
 * <ul>
 *   <li>{@link Solver#RK4}: method of lines. The age derivative is a finite-difference stencil
 *       ({@link #ageDeriv}) and time stepping is explicit RK4 for both {@code p0} and {@code L}.
 *       Explicit, so each grid step is subdivided adaptively until {@code h_sub*rate} is at most
 *       {@link #RK4_SUBSTEP_HR} (the reaction-rate stability limit); cost grows with the rates
 *       but the scheme never rejects a state on stability grounds.</li>
 *   <li>{@link Solver#SPLIT}: Strang operator splitting along characteristics. A step of length
 *       {@code h} is: aging shift by {@code h/2} ({@code f(a) <- f(a+h/2)}), reaction over
 *       {@code h} with the age-dependent rates frozen at the (now grid-aligned) characteristic
 *       midpoints, aging shift by {@code h/2}. The time step is forced to {@code 2*da} so that
 *       full-step shifts are exact index shifts ({@code deltaT} is ignored); only the partial
 *       steps at branch ends and epoch boundaries use positivity-preserving linear
 *       interpolation. The {@code L} reaction matrix is diagonal-plus-one-column
 *       ({@code -diag(alpha) + u e0^T}, the column being the age-0 coupling), so it has an O(Na)
 *       closed form given {@code p0} at the step's midpoint time. The {@code p0} reaction is
 *       a closed-form Riccati solve for {@code p0(t,0)} followed by an exponential midpoint
 *       rule per age. Unconditionally stable and positivity/[0,1] preserving, so no
 *       substepping is needed; second order in time, exact in age. Benchmarks (2026-09-22)
 *       show it is much less accurate than RK4 at equal grid spacing (about 10-100x more
 *       expensive at equal accuracy), so RK4 remains the default. Beyond the top of the (possibly
 *       truncated) age grid the shift extrapolates with the last grid value, mirroring the RK4
 *       one-sided stencil.</li>
 * </ul>
 *
 * Root conditioning: {@code conditionAt} selects what the likelihood is conditioned on.
 * {@code ORIGIN} includes the root branch (mrca to {@code maxTime}) and conditions on survival
 * from the origin. {@code MRCA} excludes the root branch and conditions on survival of both
 * mrca daughters instead; for {@code symmetric} models both daughters are age 0 at the mrca so
 * no extra information is needed. For asymmetric models one daughter continues aging from an
 * unknown pre-mrca age that can't be marginalized without a known origin time, so
 * {@code conditionAt == MRCA} additionally requires a {@code rootAge} parameter giving that
 * age explicitly (with a user-supplied prior).
 */
public class AgeDependentBirthDeathPDEModel extends AbstractModelLikelihood implements Reportable {

    /**
     * What the likelihood is conditioned on: survival from the process {@link #ORIGIN}
     * (root branch included), or survival of both {@link #MRCA} daughters (root branch
     * excluded).
     */
    public enum ConditioningPoint { ORIGIN, MRCA }

    /** Numerical scheme: method of lines with explicit {@link #RK4}, or characteristic {@link #SPLIT}. */
    public enum Solver { RK4, SPLIT }

    private static final boolean DIAG_INF = Boolean.getBoolean("beast.abd.diag");

    // Root-conditioning survival-probability guard. p0 -> 1 makes -log(1-p0) diverge to
    // +Infinity, which is a genuine unbounded runaway of the conditioned-likelihood formula
    // (not a caching artifact) once p0 lands anywhere strictly inside (1-MAX_P0_ROOT, 1) --
    // clampUnit only rejects overshoot past 1.0 exactly, so an MCMC chain can climb this
    // divergent term indefinitely via small steps while never tripping the old p0>=1.0 check.
    // MAX_P0_ROOT hard-rejects states where survival is already absurdly improbable for any
    // real observed tree; MIN_SURVIVAL floors the log argument for everything below that, so
    // the term stays large-but-finite instead of unbounded as p0 approaches the cutoff.
    private static final double MAX_P0_ROOT = 1.0 - 1e-12;
    private static final double MIN_SURVIVAL = 1e-12;

    // Explicit-RK4 stability bound for a step of size h against reaction rate r: unstable once
    // |h*r| exceeds ~2.78 (real negative eigenvalues). Confirmed 2026-08-19: without a guard,
    // an MCMC proposal that pushes birthScale/gamma into a regime where h*rate exceeds this bound
    // can make solveL/solveP0's RK4 integration overflow into a large, deceptively "good" finite
    // log-likelihood instead of the correct (tiny) one -- a one-way trapdoor an MH sampler will
    // always accept into, since the corruption only ever inflates the value, never shrinks it.
    // Rather than rejecting such states (which carved a grid-dependent hole out of the
    // posterior), each grid step is split into ceil(h*rate / RK4_SUBSTEP_HR) equal RK4
    // substeps so that every substep satisfies h_sub*rate <= RK4_SUBSTEP_HR. The target is
    // well inside the textbook bound both for safety (coupled nonlinear system, not a scalar
    // linear ODE) and for accuracy (RK4 at h*rate ~ 2 is stable but crude). Overridable via
    // -Dbeast.abd.substepHR for benchmarking.
    private static final double RK4_SUBSTEP_HR =
            Double.parseDouble(System.getProperty("beast.abd.substepHR", "1.0"));

    private static final double EPS = 1e-12;

    private final Tree tree;
    private final boolean symmetric;
    private final ConditioningPoint conditionAt;
    private final Solver solver;
    private final Parameter rootAge;

    private final Parameter birthScale;
    private final Parameter deathScale;
    private final Parameter samplingScale;
    private final Parameter extantSamplingProb;
    private final Parameter epochTimes;
    private final AgeHazard birthHazard;
    private final AgeHazard deathHazard;
    private final boolean constBirth;
    private final boolean constDeath;
    private final boolean constSampling;
    private final int numEpochs;
    private final int numBoundaries;

    private final double maxTime;
    private final int Na;
    private final int Nt;
    private final double da;
    private final double inv2da;
    private final double inv6da;
    private final double dt;
    private final double dt05;
    // p0Grid row spacing and last row index: RK4 stores p0 at half steps (dt05, 2*Nt rows),
    // SPLIT at full steps (dt == 2*da, Nt rows).
    private final double p0Step;
    private final int p0Rows;

    // Birth rate truncation
    private final double rateZeroThreshold;
    private int jLamZero;
    private int storedJLamZero;
    private int NaTrunc;
    private int storedNaTrunc;

    // Work arrays
    private double[][] branchTopL;
    private double[][] storedBranchTopL;

    private double[][] p0Grid;
    private double[][] storedP0Grid;
    private double[] p0Curr;
    private double[] p0Next;

    private double[][] birthHaz;
    private double[][] storedBirthHaz;
    private double[][] deathHaz;
    private double[][] storedDeathHaz;
    private double[] bScale;
    private double[] storedBScale;
    private double[] dScale;
    private double[] storedDScale;
    private double[] sScale;
    private double[] storedSScale;
    private double[] epBounds;
    private double[] storedEpBounds;
    // SPLIT, symmetric only: exp(-(lam_j + mu_j + psi) * dt) per epoch, the L-reaction decay
    // factor of a full step, which is time-independent within an epoch.
    private double[][] eAlphaDt;
    private double[][] storedEAlphaDt;

    // Shared tip solveL buffers
    private final int numExternal;
    private final int[] invalidTipNums;
    private final double[] invalidTipHeights;
    private final double[] invalidTipParentHeights;

    private final int[] postOrder;

    // Likelihood
    private double logLikelihood = 0.0;
    private double storedLogLikelihood = 0.0;

    // Likelihood rescaling
    private static final double SCALE_HI = 1.0e100;
    private static final double SCALE_LO = 1.0e-50;
    private double[] nodeLogScale;
    private double[] storedNodeLogScale;

    // Compute flags
    private boolean likelihoodKnown = false;
    private boolean storedLikelihoodKnown = false;
    private boolean parametersDirty = true;
    private boolean storedParametersDirty = true;
    private boolean rateStateDirty = false;
    private boolean storedStateDirty = true;

    // Node caching
    private boolean[] nodeValid;
    private boolean[] storedNodeValid;
    private int[] modifiedNodes;
    private int modifiedNodeCount = 0;

    // Parallelization
    private final int numThreads;
    private final TaskPool taskPool;

    private final int[] nodeDepth;
    private int[][] depthBuckets;
    private int[] depthBucketSizes;
    private int maxDepth;

    // Per-worker pools
    private final double[][] LPool;
    private final double[][] LmergedPool;
    private final double[][] p0BufPool;
    private final double[][] lamCurrPool;
    private final double[][] muCurrPool;
    private final double[] psiCurrPool;
    private final double[] maxRatePool;
    private final RungeKutta[] rk4Pool;
    private final double[] logScalePool;

    public AgeDependentBirthDeathPDEModel(String name,
                                                Tree tree,
                                                Parameter birthScale,
                                                AgeHazard birthHazard,
                                                Parameter deathScale,
                                                AgeHazard deathHazard,
                                                Parameter samplingScale,
                                                Parameter extantSamplingProb,
                                                Parameter epochTimes,
                                                double maxTime,
                                                double deltaA,
                                                double deltaT,
                                                boolean symmetric,
                                                ConditioningPoint conditionAt,
                                                Parameter rootAge,
                                                double rateZeroThreshold,
                                                int numThreads,
                                                Solver solver) {
        super(name);

        this.tree = tree;
        if (tree instanceof Model) {
            addModel((Model) tree);
        }
        this.symmetric = symmetric;
        this.conditionAt = conditionAt;
        this.solver = solver;
        this.rootAge = rootAge;
        if (rootAge != null) {
            addVariable(rootAge);
            rootAge.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));
        }

        this.birthScale = birthScale;
        this.constBirth = birthScale.getDimension() == 1;
        addVariable(birthScale);
        this.deathScale = deathScale;
        this.constDeath = deathScale.getDimension() == 1;
        addVariable(deathScale);
        this.samplingScale = samplingScale;
        this.constSampling = samplingScale.getDimension() == 1;
        addVariable(samplingScale);
        this.extantSamplingProb = extantSamplingProb;
        addVariable(extantSamplingProb);
        this.epochTimes = epochTimes;
        if (epochTimes != null) {
            addVariable(epochTimes);
        }
        this.numBoundaries = (epochTimes != null) ? epochTimes.getDimension() : 0;
        this.numEpochs = numBoundaries + 1;

        this.birthHazard = birthHazard;
        this.deathHazard = deathHazard;
        addModel(birthHazard);
        if (deathHazard != birthHazard) addModel(deathHazard);

        this.maxTime = maxTime;
        if (solver == Solver.SPLIT) {
            // Strang half-steps must be exact index shifts: dt == 2*m*da. DIAG: m decouples the
            // time step from the age step (m = 1 is the production setting).
            int m = Integer.getInteger("beast.abd.splitTimeMult", 1);
            this.Nt = (int) Math.ceil(maxTime / (2.0 * m * deltaA));
            this.Na = 2 * m * this.Nt;
        } else {
            this.Na = (int) Math.ceil(maxTime / deltaA);
            int NtRaw = (int) Math.ceil(maxTime / deltaT);
            this.Nt = Math.max(this.Na, NtRaw);
        }
        this.da = maxTime / this.Na;
        this.inv2da = 1.0 / (2.0 * da);
        this.inv6da = 1.0 / (6.0 * da);
        this.dt = maxTime / this.Nt;
        this.dt05 = 0.5 * dt;
        this.p0Step = (solver == Solver.SPLIT) ? this.dt : this.dt05;
        this.p0Rows = (solver == Solver.SPLIT) ? this.Nt : 2 * this.Nt;

        this.rateZeroThreshold = rateZeroThreshold;
        this.jLamZero = Na + 1;

        int totalNodes = tree.getNodeCount();

        this.branchTopL = new double[totalNodes][Na + 1];
        this.storedBranchTopL = new double[totalNodes][Na + 1];

        this.p0Grid = new double[p0Rows + 1][Na + 1];
        this.storedP0Grid = new double[p0Rows + 1][Na + 1];
        this.p0Curr = new double[Na + 1];
        this.p0Next = new double[Na + 1];

        this.postOrder = new int[totalNodes];

        this.numExternal = tree.getExternalNodeCount();
        this.invalidTipNums = new int[numExternal];
        this.invalidTipHeights = new double[numExternal];
        this.invalidTipParentHeights = new double[numExternal];

        this.birthHaz = new double[numEpochs][Na + 1];
        this.storedBirthHaz = new double[numEpochs][Na + 1];
        this.deathHaz = new double[numEpochs][Na + 1];
        this.storedDeathHaz = new double[numEpochs][Na + 1];
        this.bScale = new double[numEpochs];
        this.storedBScale = new double[numEpochs];
        this.dScale = new double[numEpochs];
        this.storedDScale = new double[numEpochs];
        this.sScale = new double[numEpochs];
        this.storedSScale = new double[numEpochs];
        this.epBounds = new double[numBoundaries];
        this.storedEpBounds = new double[numBoundaries];
        if (solver == Solver.SPLIT && symmetric) {
            this.eAlphaDt = new double[numEpochs][Na + 1];
            this.storedEAlphaDt = new double[numEpochs][Na + 1];
        }

        this.nodeLogScale = new double[totalNodes];
        this.storedNodeLogScale = new double[totalNodes];

        this.nodeValid = new boolean[totalNodes];
        this.storedNodeValid = new boolean[totalNodes];
        this.modifiedNodes = new int[totalNodes];

        this.numThreads = Math.max(1, numThreads);
        this.taskPool = (this.numThreads > 1) ? new TaskPool(this.numThreads, this.numThreads) : null;

        if (this.taskPool != null) {
            this.nodeDepth = new int[totalNodes];
            this.depthBuckets = new int[16][];
            this.depthBucketSizes = new int[16];
        } else {
            this.nodeDepth = null;
            this.depthBuckets = null;
            this.depthBucketSizes = null;
        }

        this.LPool   = new double[this.numThreads][Na + 1];
        this.LmergedPool = new double[this.numThreads][Na + 1];
        this.p0BufPool   = new double[this.numThreads][Na + 1];
        this.lamCurrPool = new double[this.numThreads][Na + 1];
        this.muCurrPool  = new double[this.numThreads][Na + 1];
        this.psiCurrPool = new double[this.numThreads];
        this.maxRatePool = new double[this.numThreads];
        this.rk4Pool     = new RungeKutta[this.numThreads];
        for (int w = 0; w < this.numThreads; w++) {
            this.rk4Pool[w] = new RungeKutta(Na + 1);
        }
        this.logScalePool = new double[this.numThreads];
    }

    /**
     * Cache scales and age-hazards
     */
    private void refreshRates() {
        for (int k = 0; k < numEpochs; k++) {
            bScale[k] = birthScale.getParameterValue(constBirth ? 0 : k);
            dScale[k] = deathScale.getParameterValue(constDeath ? 0 : k);
            sScale[k] = samplingScale.getParameterValue(constSampling ? 0 : k);
        }
        for (int k = 0; k < numBoundaries; k++) {
            epBounds[k] = epochTimes.getParameterValue(k);
        }

        for (int k = 0; k < numEpochs; k++) {
            birthHazard.evaluate(Na, da, birthHaz[k], k);
            deathHazard.evaluate(Na, da, deathHaz[k], k);
        }

        jLamZero = 0;
        for (int j = Na; j >= 0; j--) {
            double maxLam = 0.0;
            for (int k = 0; k < numEpochs; k++) {
                double v = bScale[k] * birthHaz[k][j];
                if (v > maxLam) maxLam = v;
            }
            if (maxLam >= rateZeroThreshold) {
                jLamZero = j + 1;
                break;
            }
        }

        NaTrunc = Math.min(Na, Math.max(10, Math.min(jLamZero, Na)));

        if (eAlphaDt != null) {
            for (int k = 0; k < numEpochs; k++) {
                double[] e = eAlphaDt[k];
                for (int j = 0; j <= NaTrunc; j++) {
                    e[j] = Math.exp(-(bScale[k] * birthHaz[k][j] + dScale[k] * deathHaz[k][j] + sScale[k]) * dt);
                }
            }
        }
    }

    /*
     * Current epoch rates
     */
    private void currentRates(int worker, int epoch) {
        double[] lamCurr = lamCurrPool[worker];
        double[] muCurr = muCurrPool[worker];
        double bs = bScale[epoch];
        double ds = dScale[epoch];
        double[] bHaz = birthHaz[epoch];
        double[] dHaz = deathHaz[epoch];
        double maxLam = 0.0, maxMu = 0.0;
        for (int j = 0; j <= NaTrunc; j++) {
            double lam = bs * bHaz[j];
            double mu = ds * dHaz[j];
            lamCurr[j] = lam;
            muCurr[j] = mu;
            if (lam > maxLam) maxLam = lam;
            if (mu > maxMu) maxMu = mu;
        }
        double psi = sScale[epoch];
        psiCurrPool[worker] = psi;
        // Worst-case reaction rate this epoch, used by the RK4 solver to choose the number of
        // substeps per grid step -- see RK4_SUBSTEP_HR.
        maxRatePool[worker] = maxLam + maxMu + psi;
    }

    /*
     * Age-derivative stencil for f over [0, NaTrunc], written to result.
     */
    private void ageDeriv(double[] f, double[] result) {
        final int jHi = NaTrunc - 2;
        result[0] = (-3.0 * f[0] + 4.0 * f[1] - f[2]) * inv2da;
        for (int j = 1; j <= jHi; j++) {
            result[j] = (-2.0 * f[j - 1] - 3.0 * f[j] + 6.0 * f[j + 1] - f[j + 2]) * inv6da;
        }
        // Second-order one-sided stencils, mirroring the left boundary above: the general
        // interior formula needs f[j+2], which doesn't exist for these last two points, and
        // simply dropping that term (as before) leaves stencil coefficients that don't sum to
        // zero — injecting a spurious nonzero derivative even for a perfectly flat array.
        int j1 = NaTrunc - 1;
        result[j1] = (f[j1 + 1] - f[j1 - 1]) * inv2da;
        int j2 = NaTrunc;
        result[j2] = (3.0 * f[j2] - 4.0 * f[j2 - 1] + f[j2 - 2]) * inv2da;
    }

    //===============
    // p0 solve
    //===============

    private void p0Rhs(int worker, double t, double[] p0, double[] dp0dt) {
        final double p0_0 = p0[0];
        ageDeriv(p0, dp0dt);

        double[] lamCurr = lamCurrPool[worker];
        double[] muCurr = muCurrPool[worker];
        final double psi = psiCurrPool[worker];

        if (symmetric) {
            final double p0_0sq = p0_0 * p0_0;
            for (int j = 0; j <= NaTrunc; j++) {
                double lam = lamCurr[j];
                double mu  = muCurr[j];
                dp0dt[j] += mu + lam * p0_0sq - (lam + mu + psi) * p0[j];
            }
        } else {
            for (int j = 0; j <= NaTrunc; j++) {
                double lam = lamCurr[j];
                double mu  = muCurr[j];
                double pj  = p0[j];
                dp0dt[j] += mu + lam * pj * p0_0 - (lam + mu + psi) * pj;
            }
        }
    }

    /*
     * Forward-solve p0 with half sized steps dt05. Boundary p0(0, a) = 1 - rho.
     */
    private boolean solveP0() {
        if (solver == Solver.SPLIT) return solveP0Split();
        final double rho = extantSamplingProb.getParameterValue(0);
        final double p0Init = 1.0 - rho;
        Arrays.fill(p0Curr, 0, NaTrunc + 1, p0Init);
        Arrays.fill(p0Grid[0], 0, NaTrunc + 1, p0Init);

        final int worker = 0;
        RungeKutta rk4 = rk4Pool[worker];
        RungeKutta.RhsFunction rhs = (t, y, dydt) -> p0Rhs(worker, t, y, dydt);
        int epoch = 0;
        currentRates(worker, epoch);
        double t = 0.0;

        for (int i = 1; i <= 2 * Nt; i++) {
            double tNext = i * dt05;

            if (epoch < numBoundaries && epBounds[epoch] < tNext) {
                double boundary = epBounds[epoch];
                if (boundary - t > EPS) {
                    stepP0RK4(worker, rk4, rhs, t, boundary - t);
                    t = boundary;
                }
                epoch++;
                currentRates(worker, epoch);
            }

            stepP0RK4(worker, rk4, rhs, t, tNext - t);
            t = tNext;
            System.arraycopy(p0Curr, 0, p0Grid[i], 0, NaTrunc + 1);
        }
        return true;
    }

    /* Number of equal RK4 substeps needed for a step of size h against the worker's max rate. */
    private int numSubsteps(int worker, double h) {
        double hr = h * maxRatePool[worker];
        return (hr <= RK4_SUBSTEP_HR) ? 1 : (int) Math.ceil(hr / RK4_SUBSTEP_HR);
    }

    /* Advance p0Curr by h with adaptive RK4 substepping (swaps p0Curr/p0Next). */
    private void stepP0RK4(int worker, RungeKutta rk4, RungeKutta.RhsFunction rhs, double t, double h) {
        int n = numSubsteps(worker, h);
        double hs = h / n;
        for (int k = 0; k < n; k++) {
            rk4.step(t + k * hs, hs, p0Curr, p0Next, NaTrunc + 1, rhs);
            double[] sw = p0Curr; p0Curr = p0Next; p0Next = sw;
            clampUnit(p0Curr);
        }
    }

    //===============
    // Split solver
    //===============

    /*
     * Exact aging step along characteristics, in place: x[j] <- x(a_j + h). Reads only indices
     * >= j, so ascending in-place evaluation is safe. Full steps (h a multiple of da) are pure
     * index shifts; partial steps use positivity-preserving linear interpolation. Beyond the top
     * of the grid the last value is held (constant extrapolation).
     */
    private void shiftAge(double[] x, double h) {
        final int top = NaTrunc;
        double s = h / da;
        int m = (int) Math.floor(s + 1e-9);
        double f = s - m;
        if (f < 1e-9) f = 0.0;
        if (m == 0 && f == 0.0) return;
        final double g = 1.0 - f;
        for (int j = 0; j <= top; j++) {
            int k = j + m;
            if (k >= top) {
                x[j] = x[top];
            } else if (f == 0.0) {
                x[j] = x[k];
            } else {
                x[j] = g * x[k] + f * x[k + 1];
            }
        }
    }

    /*
     * Frozen-coefficient reaction step for L over duration h with p0 evaluated at p0Mid:
     *     dL_j/dt = -alpha_j L_j + u_j L_0
     * where the coupling runs only through L_0, so L_0(s) = L_0 exp(-kappa0 s) with
     * kappa0 = alpha_0 - u_0, and each L_j has the closed form
     *     L_j(h) = e^{-alpha_j h} L_j + u_j L_0 (e^{-kappa0 h} - e^{-alpha_j h}) / (alpha_j - kappa0).
     * The quotient is evaluated via expm1 near a vanishing denominator (limit h e^{-alpha_j h}).
     * Nonnegativity is preserved: every term is a product of nonnegative factors.
     */
    private void LReactSplit(int worker, double[] L, double[] p0Mid, double h, double[] eACache) {
        double[] lamCurr = lamCurrPool[worker];
        double[] muCurr = muCurrPool[worker];
        final double psi = psiCurrPool[worker];
        final double p0_0 = p0Mid[0];
        final double x0 = L[0];

        final double lam0 = lamCurr[0];
        final double alpha0, u0;
        if (symmetric) {
            alpha0 = lam0 + muCurr[0] + psi;
            u0 = 2.0 * lam0 * p0_0;
        } else {
            alpha0 = lam0 + muCurr[0] + psi - lam0 * p0_0;
            u0 = lam0 * p0_0;
        }
        final double kappa0 = alpha0 - u0;
        final double eK = Math.exp(-kappa0 * h);

        for (int j = 0; j <= NaTrunc; j++) {
            double lam = lamCurr[j];
            double alpha, u;
            if (symmetric) {
                alpha = lam + muCurr[j] + psi;
                u = 2.0 * lam * p0_0;
            } else {
                alpha = lam + muCurr[j] + psi - lam * p0_0;
                u = lam * p0Mid[j];
            }
            double eA = (eACache != null) ? eACache[j] : Math.exp(-alpha * h);
            double d = alpha - kappa0;
            double dh = d * h;
            double q;   // (e^{-kappa0 h} - e^{-alpha h}) / (alpha - kappa0), >= 0
            if (Math.abs(dh) < 1.0) {
                q = (dh == 0.0) ? h * eA : eA * Math.expm1(dh) / d;
            } else {
                q = (eK - eA) / d;
            }
            L[j] = eA * L[j] + u * x0 * q;
        }
    }

    /*
     * Exact solution at time s of the constant-coefficient Riccati equation governing the
     * age-0 survival component during a reaction step,
     *     dp/ds = mu0 + lam0 p^2 - r0 p,   r0 = lam0 + mu0 + psi,   p(0) = pBar.
     * Roots pm <= 1 <= pp of the quadratic (real since disc >= (lam0 - mu0)^2), delta = sqrt(disc):
     *     (p - pm)/(p - pp) = C e^{-delta s},  C = (pBar - pm)/(pBar - pp).
     */
    private static double riccatiP0(double pBar, double lam0, double mu0, double r0, double s) {
        if (lam0 <= 0.0) {
            // Linear: dp/ds = mu0 - r0 p
            if (r0 <= 0.0) return pBar;
            double e = Math.exp(-r0 * s);
            return pBar * e + (mu0 / r0) * (1.0 - e);
        }
        double disc = r0 * r0 - 4.0 * lam0 * mu0;
        if (disc < 0.0) disc = 0.0;
        double delta = Math.sqrt(disc);
        double pm = (r0 - delta) / (2.0 * lam0);
        if (delta * s < 1e-8) {
            // Double root: dp/ds = lam0 (p - pm)^2
            double d = pBar - pm;
            return pm + d / (1.0 - lam0 * d * s);
        }
        double pp = (r0 + delta) / (2.0 * lam0);
        double c = (pBar - pm) / (pBar - pp);
        double ce = c * Math.exp(-delta * s);
        return (pm - pp * ce) / (1.0 - ce);
    }

    /*
     * Reaction step for p0 over duration h (between the two Strang half-shifts). The age-0
     * component obeys a closed-form Riccati equation; every other age then obeys a scalar
     * linear ODE dp_j/ds = A_j(s) - B_j(s) p_j with A_j, B_j >= 0 depending on s only through
     * p0(s,0), integrated with the exponential midpoint rule
     *     p_j(h) = e^{-B h} p_j + (A/B)(1 - e^{-B h}),   A, B at s = h/2,
     * which is second order, unconditionally stable and maps [0,1] into [0,1] since A <= B.
     */
    private void p0ReactSplit(int worker, double[] p0, double h) {
        double[] lamCurr = lamCurrPool[worker];
        double[] muCurr = muCurrPool[worker];
        final double psi = psiCurrPool[worker];

        final double lam0 = lamCurr[0];
        final double mu0 = muCurr[0];
        final double r0 = lam0 + mu0 + psi;
        final double pBar = p0[0];
        final double x0Mid = riccatiP0(pBar, lam0, mu0, r0, 0.5 * h);
        p0[0] = riccatiP0(pBar, lam0, mu0, r0, h);

        if (symmetric) {
            final double x0sq = x0Mid * x0Mid;
            for (int j = 1; j <= NaTrunc; j++) {
                double lam = lamCurr[j];
                double mu = muCurr[j];
                double A = mu + lam * x0sq;
                double B = lam + mu + psi;
                p0[j] = expMidpoint(p0[j], A, B, h);
            }
        } else {
            final double surv0 = 1.0 - x0Mid;
            for (int j = 1; j <= NaTrunc; j++) {
                double lam = lamCurr[j];
                double mu = muCurr[j];
                double B = mu + psi + lam * surv0;
                p0[j] = expMidpoint(p0[j], mu, B, h);
            }
        }
        clampUnit(p0);
    }

    /* Solution at h of dp/ds = A - B p from p, with constant A, B >= 0. */
    private static double expMidpoint(double p, double A, double B, double h) {
        double bh = B * h;
        if (bh < 1e-8) return p + (A - B * p) * h;
        double e = Math.exp(-bh);
        return p * e + (A / B) * (1.0 - e);
    }

    /*
     * Forward-solve p0 on the full-step grid (spacing dt = 2*da) with the split scheme. Boundary
     * p0(0, a) = 1 - rho. Never fails: the scheme has no stability restriction.
     */
    private boolean solveP0Split() {
        final double rho = extantSamplingProb.getParameterValue(0);
        final double p0Init = 1.0 - rho;
        Arrays.fill(p0Curr, 0, NaTrunc + 1, p0Init);
        Arrays.fill(p0Grid[0], 0, NaTrunc + 1, p0Init);

        final int worker = 0;
        int epoch = 0;
        currentRates(worker, epoch);
        double t = 0.0;

        for (int i = 1; i <= p0Rows; i++) {
            double tNext = i * p0Step;

            if (epoch < numBoundaries && epBounds[epoch] < tNext) {
                double boundary = epBounds[epoch];
                if (boundary - t > EPS) {
                    p0StepSplit(worker, p0Curr, boundary - t);
                    t = boundary;
                }
                epoch++;
                currentRates(worker, epoch);
            }

            if (tNext - t > EPS) {
                p0StepSplit(worker, p0Curr, tNext - t);
                t = tNext;
            }
            System.arraycopy(p0Curr, 0, p0Grid[i], 0, NaTrunc + 1);
        }
        return true;
    }

    /* One Strang step for p0 over duration h: half shift, reaction, half shift. */
    private void p0StepSplit(int worker, double[] p0, double h) {
        shiftAge(p0, 0.5 * h);
        p0ReactSplit(worker, p0, h);
        shiftAge(p0, 0.5 * h);
    }

    /*
     * One Strang step for L from height t to t + h: half shift, reaction with rates of the
     * current epoch and p0 at the step's midpoint time, half shift.
     */
    private void LStepSplit(int worker, double[] L, double t, double h, int epoch) {
        double[] p0Buf = p0BufPool[worker];
        double[] eACache = (eAlphaDt != null && Math.abs(h - dt) < EPS) ? eAlphaDt[epoch] : null;
        shiftAge(L, 0.5 * h);
        p0Inter(t + 0.5 * h, p0Buf);
        LReactSplit(worker, L, p0Buf, h, eACache);
        shiftAge(L, 0.5 * h);
    }

    /*
     * Split-scheme counterpart of solveL: same stepping structure (full grid steps, partial
     * steps at the branch ends and at epoch boundaries), no stability guard. Adjacent
     * half-shifts of consecutive full steps are not merged; the shift is cheap relative to the
     * exp-heavy reaction.
     */
    private void solveLSplit(int worker, double[] L, double startTime, double endTime) {
        double t = startTime;
        int epoch = getEpochIndex(t);
        currentRates(worker, epoch);

        int firstIdx = (int) Math.floor(t / dt) + 1;
        int lastIdx = (int) Math.floor(endTime / dt);

        for (int i = firstIdx; i <= lastIdx + 1; i++) {
            double tNext = (i <= lastIdx) ? i * dt : endTime;

            if (epoch < numBoundaries && epBounds[epoch] < tNext) {
                double boundary = epBounds[epoch];
                if (boundary - t > EPS) {
                    LStepSplit(worker, L, t, boundary - t, epoch);
                    if (!rescaleL(worker, L)) { L[0] = Double.NaN; return; }
                    t = boundary;
                }
                epoch++;
                currentRates(worker, epoch);
            }

            if (tNext - t > EPS) {
                LStepSplit(worker, L, t, tNext - t, epoch);
                if (!rescaleL(worker, L)) { L[0] = Double.NaN; return; }
                t = tNext;
            }
        }
    }

    /*
     * Catmull-Rom interpolation of p0 with linear fallback near the boundaries.
     */
    private void p0Inter(double t, double[] result) {
        double idx = t / p0Step;
        int lo = (int) Math.floor(idx);
        int last = p0Rows;
        if (lo < 0) lo = 0;
        if (lo >= last) {
            System.arraycopy(p0Grid[last], 0, result, 0, NaTrunc + 1);
            return;
        }
        double frac = idx - lo;
        if (lo >= 1 && lo + 2 <= last) {
            double frac2 = frac * frac;
            double frac3 = frac2 * frac;
            double wm1 = -0.5 * frac + frac2 - 0.5 * frac3;
            double w0  = 1.0 - 2.5 * frac2 + 1.5 * frac3;
            double w1  = 0.5 * frac + 2.0 * frac2 - 1.5 * frac3;
            double w2  = -0.5 * frac2 + 0.5 * frac3;
            double[] rowM1 = p0Grid[lo - 1];
            double[] row0  = p0Grid[lo];
            double[] row1  = p0Grid[lo + 1];
            double[] row2  = p0Grid[lo + 2];
            for (int j = 0; j <= NaTrunc; j++) {
                result[j] = wm1 * rowM1[j] + w0 * row0[j] + w1 * row1[j] + w2 * row2[j];
            }
        } else {
            double w0 = 1.0 - frac;
            double[] row0 = p0Grid[lo];
            double[] row1 = p0Grid[lo + 1];
            for (int j = 0; j <= NaTrunc; j++) {
                result[j] = w0 * row0[j] + frac * row1[j];
            }
        }
        // The Catmull-Rom branch above is a cubic spline through clamped [0,1] grid points --
        // the curve BETWEEN points is not guaranteed to stay in [0,1] (classic spline overshoot
        // near a sharply-curving p0(t) in extreme-parameter regions). An interpolated p0>1 fed
        // straight into LRhs's `p0AtT` (a nonphysical "probability" > 1) inflates the L-branch
        // ODE's growth term, which is the actual source of the MCMC runaway toward absurdly
        // large finite log-likelihoods (see root-conditioning MAX_P0_ROOT/MIN_SURVIVAL guard
        // above, which only catches the root term and was insufficient on its own). Clamp here
        // so every consumer of p0Inter sees a valid probability.
        for (int j = 0; j <= NaTrunc; j++) {
            result[j] = Math.min(Math.max(result[j], 0.0), 1.0);
        }
    }

    //===============
    // L solve
    //===============

    private void LRhs(int worker, double[] L, double[] p0AtT, double[] dLdt) {
        final double L_0 = L[0];
        final double p0_0 = p0AtT[0];
        ageDeriv(L, dLdt);

        double[] lamCurr = lamCurrPool[worker];
        double[] muCurr = muCurrPool[worker];
        final double psi = psiCurrPool[worker];

        if (symmetric) {
            final double c = 2.0 * p0_0 * L_0;
            for (int j = 0; j <= NaTrunc; j++) {
                double lam = lamCurr[j];
                double r   = lam + muCurr[j] + psi;
                dLdt[j] += c * lam - r * L[j];
            }
        } else {
            for (int j = 0; j <= NaTrunc; j++) {
                double lam = lamCurr[j];
                double r   = lam + muCurr[j] + psi;
                dLdt[j] += lam * (L[j] * p0_0 + p0AtT[j] * L_0) - r * L[j];
            }
        }
    }

    /*
     * Solve L over [startTime, endTime]
     */
    private static final boolean DIAG3 = Boolean.getBoolean("beast.abd.diag3");

    private String argmaxStr(double[] L) {
        int argmax = 0;
        double max = L[0];
        for (int j = 1; j <= NaTrunc; j++) {
            if (L[j] > max) { max = L[j]; argmax = j; }
        }
        return "argmaxJ=" + argmax + " maxVal=" + max + " maxAge=" + (argmax * da);
    }

    private void solveL(int worker, double[] L, double startTime, double endTime) {
        if (endTime - startTime < EPS) return;
        if (solver == Solver.SPLIT) { solveLSplit(worker, L, startTime, endTime); return; }
        if (DIAG3) System.err.println("solveL-START worker=" + worker + " startTime=" + startTime
                + " endTime=" + endTime + " L0=" + L[0] + " " + argmaxStr(L));

        double[] p0Buf = p0BufPool[worker];
        RungeKutta rk4 = rk4Pool[worker];

        RungeKutta.RhsFunction rhsInterp = (t, y, dydt) -> {
            p0Inter(t, p0Buf);
            LRhs(worker, y, p0Buf, dydt);
        };

        RungeKutta.RhsFunction rhsGrid = (t, y, dydt) -> {
            int idx = (int) Math.round(t / p0Step);
            LRhs(worker, y, p0Grid[idx], dydt);
        };

        double t = startTime;
        int epoch = getEpochIndex(t);
        currentRates(worker, epoch);

        int firstIdx = (int) Math.floor(t / dt) + 1;
        int lastIdx = (int) Math.floor(endTime / dt);
        boolean tOnGrid = Math.abs(t - (firstIdx - 1) * dt) < EPS;

        for (int i = firstIdx; i <= lastIdx + 1; i++) {
            double tNext = (i <= lastIdx) ? i * dt : endTime;

            if (epoch < numBoundaries && epBounds[epoch] < tNext) {
                double boundary = epBounds[epoch];
                if (boundary - t > EPS) {
                    stepLRK4(worker, rk4, rhsInterp, t, boundary - t, L);
                    clampNonNeg(L);
                    if (!rescaleL(worker, L)) { L[0] = Double.NaN; return; }
                    t = boundary;
                    tOnGrid = false;
                    if (DIAG3) System.err.println("solveL-step(bound) t=" + t + " L0=" + L[0]
                            + " logScale=" + logScalePool[worker] + " " + argmaxStr(L));
                }
                epoch++;
                currentRates(worker, epoch);
            }

            if (tNext - t > EPS) {
                boolean nextOnGrid = (i <= lastIdx);
                boolean stepOnGrid = tOnGrid && nextOnGrid && Math.abs((tNext - t) - dt) < EPS;
                int n = numSubsteps(worker, tNext - t);
                if (n == 1) {
                    rk4.step(t, tNext - t, L, L, NaTrunc + 1, stepOnGrid ? rhsGrid : rhsInterp);
                } else {
                    // Substeps evaluate p0 off the half-step grid, so always interpolate.
                    stepLRK4(worker, rk4, rhsInterp, t, tNext - t, L);
                }
                clampNonNeg(L);
                if (!rescaleL(worker, L)) { L[0] = Double.NaN; return; }
                t = tNext;
                tOnGrid = nextOnGrid;
                if (DIAG3) System.err.println("solveL-step t=" + t + " L0=" + L[0]
                        + " logScale=" + logScalePool[worker] + " " + argmaxStr(L));
            }
        }
    }

    /* Advance L in place by h with adaptive RK4 substepping. */
    private void stepLRK4(int worker, RungeKutta rk4, RungeKutta.RhsFunction rhs, double t, double h, double[] L) {
        int n = numSubsteps(worker, h);
        double hs = h / n;
        for (int k = 0; k < n; k++) {
            rk4.step(t + k * hs, hs, L, L, NaTrunc + 1, rhs);
        }
    }

    private boolean rescaleL(int worker, double[] L) {
        double maxL = 0.0;
        for (int j = 0; j <= NaTrunc; j++) {
            double v = L[j];
            if (!Double.isFinite(v)) return false;
            if (v > maxL) maxL = v;
        }
        if (maxL == 0.0) return false;
        if (maxL > SCALE_HI || maxL < SCALE_LO) {
            double inv = 1.0 / maxL;
            for (int j = 0; j <= NaTrunc; j++) {
                L[j] *= inv;
            }
            logScalePool[worker] += Math.log(maxL);
        }
        return true;
    }

    /*
     * Compute L over all external branches in one pass. Tips that share an identical
     * (tipHeight, parentHeight) pair are true cherries: their L integral covers the exact same
     * range, so the solved buffer is reused verbatim. Tips are NOT chained across different
     * parentHeights (even within the same tipHeight group): the RK4 grid stepping in solveL
     * switches to interpolated rates at a non-grid-aligned checkpoint, so splitting one
     * grid-aligned step into two checkpoint-bounded steps is not numerically equivalent to
     * taking it whole, even though the underlying continuous L(t) trajectory is tip-independent.
     * Reusing a buffer across unrelated tips therefore made a tip's cached value depend on
     * which other tips happened to be invalidated alongside it in that round.
     */
    private void solveLTips() {
        int numInvalid = 0;

        for (int i = 0; i < numExternal; i++) {
            NodeRef node = tree.getExternalNode(i);
            int nodeNum = node.getNumber();
            if (nodeValid[nodeNum]) continue;

            double tipHeight = tree.getNodeHeight(node);
            double parentHeight = tree.getNodeHeight(tree.getParent(node));

            invalidTipNums[numInvalid] = nodeNum;
            invalidTipHeights[numInvalid] = tipHeight;
            invalidTipParentHeights[numInvalid] = parentHeight;
            numInvalid++;
        }

        if (numInvalid == 0) return;

        sortByTipThenParent(invalidTipNums, invalidTipHeights, invalidTipParentHeights, numInvalid);

        final int worker = 0;
        double[] L = LPool[worker];

        double groupTipHeight = Double.NaN;
        double groupParentHeight = Double.NaN;

        for (int i = 0; i < numInvalid; i++) {
            double tipHeight = invalidTipHeights[i];
            double parentHeight = invalidTipParentHeights[i];
            int tipNum = invalidTipNums[i];

            if (tipHeight != groupTipHeight || parentHeight != groupParentHeight) {
                Arrays.fill(L, 0, NaTrunc + 1, 1.0);
                logScalePool[worker] = 0.0;
                solveL(worker, L, tipHeight, parentHeight);
                groupTipHeight = tipHeight;
                groupParentHeight = parentHeight;
            }

            nodeLogScale[tipNum] = logScalePool[worker];
            System.arraycopy(L, 0, branchTopL[tipNum], 0, NaTrunc + 1);
            modifiedNodes[modifiedNodeCount++] = tipNum;
        }
    }

    /*
     * Solve L for a given internal branch
     */
    private void solveLInternal(int worker, int rootNum) {
        NodeRef root = tree.getNode(rootNum);
        double rootHeight = tree.getNodeHeight(root);
        double endTime = tree.isRoot(root) ? maxTime : tree.getNodeHeight(tree.getParent(root));

        int left = tree.getChild(root, 0).getNumber();
        int right = tree.getChild(root, 1).getNumber();

        double leftL_0 = branchTopL[left][0];
        double rightL_0 = branchTopL[right][0];

        currentRates(worker, getEpochIndex(rootHeight));
        double[] Lmerged = LmergedPool[worker];
        double[] lamCurr = lamCurrPool[worker];
        for (int j = 0; j <= NaTrunc; j++) {
            double lam = lamCurr[j];
            if (symmetric) {
                Lmerged[j] = lam * leftL_0 * rightL_0;
            } else {
                Lmerged[j] = lam * (branchTopL[left][j] * rightL_0
                        + leftL_0 * branchTopL[right][j]);
            }
        }

        logScalePool[worker] = 0.0;
        double[] dst = branchTopL[rootNum];
        System.arraycopy(Lmerged, 0, dst, 0, NaTrunc + 1);
        solveL(worker, dst, rootHeight, endTime);
        nodeLogScale[rootNum] = logScalePool[worker];
    }

    // =========================
    // Likelihood computation
    // =========================
    private double calculateLogLikelihood() {
        double rootH = tree.getNodeHeight(tree.getRoot());

        if (maxTime <= rootH) {
            if (DIAG_INF) System.err.println("ABD -Inf [rootHeight]: rootHeight="
                    + rootH + " >= maxTime=" + maxTime);
            return Double.NEGATIVE_INFINITY;
        }

        if (parametersDirty) {
            refreshRates();
            boolean p0Ok = solveP0();
            rateStateDirty = true;
            parametersDirty = false;
            if (!p0Ok) {
                if (DIAG_INF) System.err.println("ABD -Inf [p0Solve]: p0 solve failed");
                return Double.NEGATIVE_INFINITY;
            }
        }

        modifiedNodeCount = 0;

        solveLTips();

        TreeUtils.postOrderTraversalList(tree, postOrder);

        if (taskPool == null) {
            for (int nodeNum : postOrder) {
                NodeRef node = tree.getNode(nodeNum);
                if (nodeValid[nodeNum] || tree.isExternal(node)
                        || (conditionAt == ConditioningPoint.MRCA && tree.isRoot(node))) continue;
                solveLInternal(0, nodeNum);
                modifiedNodes[modifiedNodeCount++] = nodeNum;
            }
        } else {
            buildDepthBuckets();
            AtomicInteger modifiedSlot = new AtomicInteger(modifiedNodeCount);
            for (int lvl = 1; lvl <= maxDepth; lvl++) {
                int n = depthBucketSizes[lvl];
                if (n == 0) continue;
                int[] bucket = depthBuckets[lvl];

                if (n == 1) {
                    solveLInternal(0, bucket[0]);
                    modifiedNodes[modifiedSlot.getAndIncrement()] = bucket[0];
                    continue;
                }

                final int bucketSize = n;
                final int[] bucketRef = bucket;
                final AtomicInteger nextItem = new AtomicInteger(0);
                taskPool.fork((task, thread) -> {
                    int idx;
                    while ((idx = nextItem.getAndIncrement()) < bucketSize) {
                        solveLInternal(thread, bucketRef[idx]);
                        modifiedNodes[modifiedSlot.getAndIncrement()] = bucketRef[idx];
                    }
                });
            }
            modifiedNodeCount = modifiedSlot.get();
        }

        Arrays.fill(nodeValid, true);

        double logLik = 0.0;
        double totalLogScale = 0.0;
        for (double v : nodeLogScale) totalLogScale += v;

        // Per-tip sampling factor: psi(t_i) for serial, rho for extant
        double logSamplingFactor = 0.0;
        double rho = extantSamplingProb.getParameterValue(0);
        double logRho = (rho > 0.0) ? Math.log(rho) : Double.NEGATIVE_INFINITY;
        for (int i = 0; i < numExternal; i++) {
            NodeRef tip = tree.getExternalNode(i);
            double th = tree.getNodeHeight(tip);
            if (th > 0.0) {
                double psi = sScale[getEpochIndex(th)];
                if (psi <= 0.0) {
                    if (DIAG_INF) System.err.println("ABD -Inf [psi<=0]: psi=" + psi
                            + " tipHeight=" + th);
                    return Double.NEGATIVE_INFINITY;
                }
                logSamplingFactor += Math.log(psi);
            } else {
                logSamplingFactor += logRho;
            }
        }

        NodeRef root = tree.getRoot();
        int rootNum = root.getNumber();
        if (conditionAt == ConditioningPoint.MRCA) {
            double rootHeight = tree.getNodeHeight(root);
            int left = tree.getChild(root, 0).getNumber();
            int right = tree.getChild(root, 1).getNumber();

            double leftL_0 = branchTopL[left][0];
            double rightL_0 = branchTopL[right][0];
            if (!Double.isFinite(leftL_0) || !Double.isFinite(rightL_0)
                    || leftL_0 <= 0.0 || rightL_0 <= 0.0) {
                if (DIAG_INF) System.err.println("ABD -Inf [rootChildL]: leftL_0="
                        + leftL_0 + " rightL_0=" + rightL_0 + " rootHeight=" + rootHeight);
                return Double.NEGATIVE_INFINITY;
            }

            double[] p0Buf = p0BufPool[0];
            p0Inter(rootHeight, p0Buf);

            double p0Root = p0Buf[0];
            if (!Double.isFinite(p0Root) || p0Root >= MAX_P0_ROOT) {
                if (DIAG_INF) System.err.println("ABD -Inf [p0Root]: p0Root="
                        + p0Root + " rootHeight=" + rootHeight);
                return Double.NEGATIVE_INFINITY;
            }

            double logLPartial;
            double logSurvival;
            if (symmetric) {
                logLPartial = Math.log(leftL_0) + Math.log(rightL_0);
                logSurvival = -2.0 * Math.log(Math.max(1.0 - p0Root, MIN_SURVIVAL));
            } else {
                // Asymmetric mrca-conditioning: one daughter continues aging from an unknown
                // pre-mrca age, given here via rootAge (with its own prior) since it cannot be
                // marginalized without a known origin time. Symmetrized over which daughter is
                // the continuing one, mirroring solveLInternal's Lmerged but without the
                // birth-rate factor (the split is conditioned on as given, not scored).
                double aStar = rootAge.getParameterValue(0);
                if (!Double.isFinite(aStar) || aStar < 0.0 || aStar > NaTrunc * da) {
                    if (DIAG_INF) System.err.println("ABD -Inf [rootAge]: rootAge="
                            + aStar + " maxAge=" + (NaTrunc * da));
                    return Double.NEGATIVE_INFINITY;
                }

                double leftL_a = interpAge(branchTopL[left], aStar);
                double rightL_a = interpAge(branchTopL[right], aStar);
                double partial = leftL_0 * rightL_a + leftL_a * rightL_0;
                if (!Double.isFinite(partial) || partial <= 0.0) {
                    if (DIAG_INF) System.err.println("ABD -Inf [rootChildL_a]: leftL_a="
                            + leftL_a + " rightL_a=" + rightL_a + " rootAge=" + aStar);
                    return Double.NEGATIVE_INFINITY;
                }

                double p0RootA = interpAge(p0Buf, aStar);
                if (!Double.isFinite(p0RootA) || p0RootA >= MAX_P0_ROOT) {
                    if (DIAG_INF) System.err.println("ABD -Inf [p0RootA]: p0RootA="
                            + p0RootA + " rootAge=" + aStar);
                    return Double.NEGATIVE_INFINITY;
                }

                logLPartial = Math.log(partial);
                logSurvival = -Math.log(Math.max(1.0 - p0Root, MIN_SURVIVAL))
                        - Math.log(Math.max(1.0 - p0RootA, MIN_SURVIVAL));
                if (Boolean.getBoolean("beast.abd.diag2")) {
                    System.err.println("ABD-BREAKDOWN-ASYM: aStar=" + aStar
                            + " leftL_0=" + leftL_0 + " rightL_0=" + rightL_0
                            + " leftL_a=" + leftL_a + " rightL_a=" + rightL_a
                            + " partial=" + partial + " p0Root=" + p0Root + " p0RootA=" + p0RootA
                            + " logLPartial=" + logLPartial + " logSurvival=" + logSurvival);
                }
            }
            logLik += logLPartial + logSurvival + totalLogScale + logSamplingFactor;
            if (Boolean.getBoolean("beast.abd.diag2")) {
                System.err.println("ABD-BREAKDOWN-MRCA: logLPartial=" + logLPartial
                        + " logSurvival=" + logSurvival + " p0Root=" + p0Root
                        + " leftL_0=" + leftL_0 + " rightL_0=" + rightL_0);
            }
        } else {
            double LRoot_0 = branchTopL[rootNum][0];
            if (!Double.isFinite(LRoot_0) || LRoot_0 <= 0.0) {
                if (DIAG_INF) System.err.println("ABD -Inf [rootL]: LRoot_0=" + LRoot_0);
                return Double.NEGATIVE_INFINITY;
            }

            double p0Origin = p0Grid[p0Rows][0];
            if (!Double.isFinite(p0Origin) || p0Origin >= MAX_P0_ROOT) {
                if (DIAG_INF) System.err.println("ABD -Inf [p0Origin]: p0Origin=" + p0Origin);
                return Double.NEGATIVE_INFINITY;
            }

            double logLPartial = Math.log(branchTopL[rootNum][0]);
            double logSurvival = -1.0 * Math.log(Math.max(1.0 - p0Origin, MIN_SURVIVAL));
            logLik += logLPartial + logSurvival + totalLogScale + logSamplingFactor;
        }

        if (!Double.isFinite(logLik)) {
            if (DIAG_INF) System.err.println("ABD -Inf [nonFiniteLogLik]: logLik=" + logLik
                    + " totalLogScale=" + totalLogScale);
            return Double.NEGATIVE_INFINITY;
        }

        if (Boolean.getBoolean("beast.abd.diag2")) {
            System.err.println("ABD-BREAKDOWN: logLik=" + logLik
                    + " totalLogScale=" + totalLogScale + " logSamplingFactor=" + logSamplingFactor
                    + " NaTrunc=" + NaTrunc + " jLamZero=" + jLamZero);
        }

        return logLik;
    }

    // =======
    // Utils
    // =======

    private int getEpochIndex(double t) {
        int epoch = 0;
        while (epoch < numBoundaries && t >= epBounds[epoch]) {
            epoch++;
        }
        return epoch;
    }

    /*
     * Linear interpolation of a row indexed by age (spacing da) at age a. Clamps to
     * [0, NaTrunc*da]. Used for both branchTopL rows and p0 rows, which share the same
     * age indexing.
     */
    private double interpAge(double[] row, double a) {
        double idx = a / da;
        if (idx <= 0.0) return row[0];
        if (idx >= NaTrunc) return row[NaTrunc];
        int lo = (int) Math.floor(idx);
        double frac = idx - lo;
        return (1.0 - frac) * row[lo] + frac * row[lo + 1];
    }

    /*
     * Insertion sort of three parallel arrays by (tipKeys, parentKeys) lexicographically.
     */
    private static void sortByTipThenParent(int[] values, double[] tipKeys, double[] parentKeys, int n) {
        for (int k = 1; k < n; k++) {
            int v = values[k];
            double tk = tipKeys[k];
            double pk = parentKeys[k];
            int m = k - 1;
            while (m >= 0 && (tipKeys[m] > tk || (tipKeys[m] == tk && parentKeys[m] > pk))) {
                values[m + 1] = values[m];
                tipKeys[m + 1] = tipKeys[m];
                parentKeys[m + 1] = parentKeys[m];
                m--;
            }
            values[m + 1] = v;
            tipKeys[m + 1] = tk;
            parentKeys[m + 1] = pk;
        }
    }

    private void clampNonNeg(double[] x) {
        for (int j = 0; j <= NaTrunc; j++) {
            x[j] = Math.max(x[j], 0.0);
        }
    }

    private void clampUnit(double[] x) {
        for (int j = 0; j <= NaTrunc; j++) {
            x[j] = Math.min(Math.max(x[j], 0.0), 1.0);
        }
    }

    // ========================
    // Parallelization utils
    // ========================

    private void buildDepthBuckets() {
        Arrays.fill(depthBucketSizes, 0);
        maxDepth = 0;

        for (int nodeNum : postOrder) {
            NodeRef node = tree.getNode(nodeNum);
            if (tree.isExternal(node)) {
                nodeDepth[nodeNum] = 0;
                continue;
            }

            int leftLvl = nodeDepth[tree.getChild(node, 0).getNumber()];
            int rightLvl = nodeDepth[tree.getChild(node, 1).getNumber()];
            int lvl = 1 + Math.max(leftLvl, rightLvl);
            nodeDepth[nodeNum] = lvl;

            if (nodeValid[nodeNum]) continue;
            if (conditionAt == ConditioningPoint.MRCA && tree.isRoot(node)) continue;

            ensureDepthCapacity(lvl);
            int sz = depthBucketSizes[lvl];
            int[] bucket = depthBuckets[lvl];
            if (sz >= bucket.length) {
                bucket = Arrays.copyOf(bucket, bucket.length * 2);
                depthBuckets[lvl] = bucket;
            }
            bucket[sz] = nodeNum;
            depthBucketSizes[lvl] = sz + 1;
            if (lvl > maxDepth) maxDepth = lvl;
        }
    }

    private void ensureDepthCapacity(int level) {
        if (level < depthBuckets.length) {
            if (depthBuckets[level] == null) depthBuckets[level] = new int[16];
            return;
        }
        int newSize = Math.max(level + 1, depthBuckets.length * 2);
        depthBuckets = Arrays.copyOf(depthBuckets, newSize);
        depthBucketSizes = Arrays.copyOf(depthBucketSizes, newSize);
        depthBuckets[level] = new int[16];
    }

    // =========================
    // MCMC state management
    // =========================

    public Model getModel() {
        return this;
    }

    public double getLogLikelihood() {
        if (!likelihoodKnown) {
            logLikelihood = calculateLogLikelihood();
            likelihoodKnown = true;
        }
        return logLikelihood;
    }

    /*
     * Invalidate a node's branchTopL cache and everything whose cache could
     * legitimately depend on it.
     */
    private void invalidateNode(NodeRef node) {
        int nChildren = tree.getChildCount(node);
        for (int i = 0; i < nChildren; i++) {
            nodeValid[tree.getChild(node, i).getNumber()] = false;
        }

        NodeRef current = node;
        while (current != null) {
            nodeValid[current.getNumber()] = false;
            if (tree.isRoot(current)) break;
            NodeRef parent = tree.getParent(current);
            int siblingCount = tree.getChildCount(parent);
            for (int i = 0; i < siblingCount; i++) {
                NodeRef sibling = tree.getChild(parent, i);
                if (sibling.getNumber() != current.getNumber()) {
                    nodeValid[sibling.getNumber()] = false;
                }
            }
            current = parent;
        }
    }

    private void invalidateAllNodes() {
        Arrays.fill(nodeValid, false);
    }

    public void makeDirty() {
        likelihoodKnown = false;
        parametersDirty = true;
        invalidateAllNodes();
    }

    protected void handleModelChangedEvent(Model model, Object object, int index) {
        if (model == tree) {
            if (object instanceof TreeChangedEvent) {
                TreeChangedEvent event = (TreeChangedEvent) object;
                if (event.getNode() != null) {
                    invalidateNode(event.getNode());
                } else {
                    invalidateAllNodes();
                }
            } else {
                invalidateAllNodes();
            }
            likelihoodKnown = false;
        } else {
            makeDirty();
        }
    }

    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        makeDirty();
    }

    protected void storeState() {
        storedLogLikelihood = logLikelihood;
        storedLikelihoodKnown = likelihoodKnown;
        storedParametersDirty = parametersDirty;
        storedJLamZero = jLamZero;
        storedNaTrunc = NaTrunc;

        if (storedStateDirty || rateStateDirty) {
            int totalCols = p0Rows + 1;
            for (int i = 0; i < totalCols; i++) {
                System.arraycopy(p0Grid[i], 0, storedP0Grid[i], 0, NaTrunc + 1);
            }
            for (int k = 0; k < numEpochs; k++) {
                System.arraycopy(birthHaz[k], 0, storedBirthHaz[k], 0, Na + 1);
                System.arraycopy(deathHaz[k], 0, storedDeathHaz[k], 0, Na + 1);
                if (eAlphaDt != null) System.arraycopy(eAlphaDt[k], 0, storedEAlphaDt[k], 0, Na + 1);
            }
            System.arraycopy(bScale, 0, storedBScale, 0, bScale.length);
            System.arraycopy(dScale, 0, storedDScale, 0, dScale.length);
            System.arraycopy(sScale, 0, storedSScale, 0, sScale.length);
            System.arraycopy(epBounds, 0, storedEpBounds, 0, numBoundaries);
            rateStateDirty = false;
        }

        if (storedStateDirty) {
            for (int i = 0; i < branchTopL.length; i++) {
                System.arraycopy(branchTopL[i], 0, storedBranchTopL[i], 0, NaTrunc + 1);
            }
            System.arraycopy(nodeLogScale, 0, storedNodeLogScale, 0, nodeLogScale.length);
            storedStateDirty = false;
        } else {
            for (int i = 0; i < modifiedNodeCount; i++) {
                int n = modifiedNodes[i];
                System.arraycopy(branchTopL[n], 0, storedBranchTopL[n], 0, NaTrunc + 1);
                storedNodeLogScale[n] = nodeLogScale[n];
            }
        }
        modifiedNodeCount = 0;

        System.arraycopy(nodeValid, 0, storedNodeValid, 0, nodeValid.length);
    }

    protected void restoreState() {
        logLikelihood = storedLogLikelihood;
        likelihoodKnown = storedLikelihoodKnown;
        parametersDirty = storedParametersDirty;
        jLamZero = storedJLamZero;
        NaTrunc = storedNaTrunc;

        if (rateStateDirty) {
            double[][] tmp2D;
            double[] tmpD;
            boolean[] tmpB;

            tmp2D = p0Grid; p0Grid = storedP0Grid; storedP0Grid = tmp2D;
            tmp2D = birthHaz; birthHaz = storedBirthHaz; storedBirthHaz = tmp2D;
            tmp2D = deathHaz; deathHaz = storedDeathHaz; storedDeathHaz = tmp2D;
            tmp2D = eAlphaDt; eAlphaDt = storedEAlphaDt; storedEAlphaDt = tmp2D;
            tmpD = bScale; bScale = storedBScale; storedBScale = tmpD;
            tmpD = dScale; dScale = storedDScale; storedDScale = tmpD;
            tmpD = sScale; sScale = storedSScale; storedSScale = tmpD;
            tmpD = epBounds; epBounds = storedEpBounds; storedEpBounds = tmpD;
            tmp2D = branchTopL; branchTopL = storedBranchTopL; storedBranchTopL = tmp2D;
            tmpD = nodeLogScale; nodeLogScale = storedNodeLogScale; storedNodeLogScale = tmpD;
            tmpB = nodeValid; nodeValid = storedNodeValid; storedNodeValid = tmpB;

            storedStateDirty = true;
        } else {
            for (int i = 0; i < modifiedNodeCount; i++) {
                int n = modifiedNodes[i];
                System.arraycopy(storedBranchTopL[n], 0, branchTopL[n], 0, NaTrunc + 1);
                nodeLogScale[n] = storedNodeLogScale[n];
            }
            System.arraycopy(storedNodeValid, 0, nodeValid, 0, nodeValid.length);
        }

        rateStateDirty = false;
        modifiedNodeCount = 0;
    }

    protected void acceptState() {
    }

    public String getReport() {
        getLogLikelihood();
        return "logLikelihood: " + logLikelihood
                + " (NaTrunc=" + NaTrunc + "/" + Na
                + ", jLamZero=" + jLamZero
                + ", numThreads=" + numThreads
                + ", rho=" + extantSamplingProb.getParameterValue(0)
                + ", solver=" + solver
                + ", p0Origin=" + p0Grid[p0Rows][0] + ")\n";
    }

    public String toString() {
        return Double.toString(getLogLikelihood());
    }
}
