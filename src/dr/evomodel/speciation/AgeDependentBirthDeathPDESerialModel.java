package dr.evomodel.speciation;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeUtils;
import dr.evomodel.tree.TreeChangedEvent;
import dr.inference.model.*;
import dr.math.RungeKutta;
import dr.xml.Reportable;
import dr.util.TaskPool;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.Arrays;

/**
 * @author Frederik M. Andersen
 *
 * Serially-sampled variant of {@link AgeDependentBirthDeathPDEModel}.
 * Tips may have height > 0 (fossils / serial samples). The model adds a time-only serial sampling
 * rate psi(t) (skyline) and an extant sampling probability rho.
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
 * Conditioning expressions are unchanged in form relative to the ultrametric model; only p0
 * carries new semantics.
 *
 * Structurally a sibling of {@link AgeDependentBirthDeathPDEModel}: field names, helper names
 * and method signatures are kept identical wherever behavior is identical, so the two classes
 * can be merged in the future with a minimal diff.
 */
public class AgeDependentBirthDeathPDESerialModel extends AbstractModelLikelihood implements Reportable {
    private static final boolean DIAG_INF = Boolean.getBoolean("beast.abd.diag");

    private static final double EPS = 1e-12;

    private final Tree tree;
    private final boolean symmetric;
    private final boolean excludeRootBranch;

    private final Parameter birthScale;
    private final Parameter birthShape;
    private final Parameter deathScale;
    private final Parameter deathShape;
    private final Parameter samplingScale;
    private final Parameter extantSamplingProb;
    private final Parameter epochTimes;
    private final boolean constBirth;
    private final boolean constDeath;
    private final boolean constSampling;
    private final boolean constBirthShape;
    private final boolean constDeathShape;
    private final int numEpochs;
    private final int numBoundaries;

    private final double originTime;
    private final int Na;
    private final int Nt;
    private final double da;
    private final double inv2da;
    private final double inv6da;
    private final double dt;
    private final double dt05;

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
    private final RungeKutta[] rk4Pool;
    private final double[] logScalePool;

    public AgeDependentBirthDeathPDESerialModel(String name,
                                                Tree tree,
                                                Parameter birthScale,
                                                Parameter birthShape,
                                                Parameter deathScale,
                                                Parameter deathShape,
                                                Parameter samplingScale,
                                                Parameter extantSamplingProb,
                                                Parameter epochTimes,
                                                double originTime,
                                                int Na,
                                                int Nt,
                                                boolean symmetric,
                                                boolean excludeRootBranch,
                                                double rateZeroThreshold,
                                                int numThreads) {
        super(name);

        this.tree = tree;
        if (tree instanceof Model) {
            addModel((Model) tree);
        }
        this.symmetric = symmetric;
        this.excludeRootBranch = excludeRootBranch;

        this.birthScale = birthScale;
        this.constBirth = birthScale.getDimension() == 1;
        addVariable(birthScale);
        this.birthShape = birthShape;
        this.constBirthShape = birthShape.getDimension() == 2;
        addVariable(birthShape);
        this.deathScale = deathScale;
        this.constDeath = deathScale.getDimension() == 1;
        addVariable(deathScale);
        this.deathShape = deathShape;
        this.constDeathShape = deathShape.getDimension() == 2;
        addVariable(deathShape);
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

        this.originTime = originTime;
        this.Na = Na;
        this.Nt = Math.max(Na, Nt);
        this.da = originTime / Na;
        this.inv2da = 1.0 / (2.0 * da);
        this.inv6da = 1.0 / (6.0 * da);
        this.dt = originTime / this.Nt;
        this.dt05 = 0.5 * dt;

        // Birth rate truncation
        this.rateZeroThreshold = rateZeroThreshold;
        this.jLamZero = Na + 1;

        // Work arrays
        int totalNodes = tree.getNodeCount();
        int Nt2 = 2 * this.Nt;

        this.branchTopL = new double[totalNodes][Na + 1];
        this.storedBranchTopL = new double[totalNodes][Na + 1];

        this.p0Grid = new double[Nt2 + 1][Na + 1];
        this.storedP0Grid = new double[Nt2 + 1][Na + 1];
        this.p0Curr = new double[Na + 1];
        this.p0Next = new double[Na + 1];

        this.postOrder = new int[totalNodes];

        // Shared tip solveL buffers
        this.numExternal = tree.getExternalNodeCount();
        this.invalidTipNums = new int[numExternal];
        this.invalidTipHeights = new double[numExternal];
        this.invalidTipParentHeights = new double[numExternal];

        // Hazards
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

        // Likelihood rescaling
        this.nodeLogScale = new double[totalNodes];
        this.storedNodeLogScale = new double[totalNodes];

        // Node caching
        this.nodeValid = new boolean[totalNodes];
        this.storedNodeValid = new boolean[totalNodes];
        this.modifiedNodes = new int[totalNodes];

        // Parallelization
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

        // Per-worker pools
        this.LPool   = new double[this.numThreads][Na + 1];
        this.LmergedPool = new double[this.numThreads][Na + 1];
        this.p0BufPool   = new double[this.numThreads][Na + 1];
        this.lamCurrPool = new double[this.numThreads][Na + 1];
        this.muCurrPool  = new double[this.numThreads][Na + 1];
        this.psiCurrPool = new double[this.numThreads];
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
            int bOff = constBirthShape ? 0 : 2 * k;
            int dOff = constDeathShape ? 0 : 2 * k;
            double bShapeR     = birthShape.getParameterValue(bOff);
            double bShapeGamma = birthShape.getParameterValue(bOff + 1);
            double dShapeR     = deathShape.getParameterValue(dOff);
            double dShapeGamma = deathShape.getParameterValue(dOff + 1);

            double bShapeLin = bShapeR * bShapeGamma;
            double dShapeLin = dShapeR * dShapeGamma;
            double[] bHazK = birthHaz[k];
            double[] dHazK = deathHaz[k];
            for (int j = 0; j <= Na; j++) {
                double a = j * da;
                bHazK[j] = (1.0 + bShapeLin * a) * Math.exp(-bShapeGamma * a);
                dHazK[j] = (1.0 + dShapeLin * a) * Math.exp(-dShapeGamma * a);
            }
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
        for (int j = 0; j <= NaTrunc; j++) {
            lamCurr[j] = bs * bHaz[j];
            muCurr[j] = ds * dHaz[j];
        }
        psiCurrPool[worker] = sScale[epoch];
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
        int j1 = NaTrunc - 1;
        result[j1] = (-2.0 * f[j1 - 1] - 3.0 * f[j1] + 6.0 * f[j1 + 1]) * inv6da;
        int j2 = NaTrunc;
        result[j2] = (-2.0 * f[j2 - 1] - 3.0 * f[j2]) * inv6da;
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
    private void solveP0() {
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
                    rk4.step(t, boundary - t, p0Curr, p0Next, NaTrunc + 1, rhs);
                    double[] sw = p0Curr; p0Curr = p0Next; p0Next = sw;
                    clampUnit(p0Curr);
                    t = boundary;
                }
                epoch++;
                currentRates(worker, epoch);
            }

            rk4.step(t, tNext - t, p0Curr, p0Next, NaTrunc + 1, rhs);
            double[] sw = p0Curr; p0Curr = p0Next; p0Next = sw;
            clampUnit(p0Curr);
            t = tNext;
            System.arraycopy(p0Curr, 0, p0Grid[i], 0, NaTrunc + 1);
        }
    }

    /*
     * Catmull-Rom interpolation of p0 with linear fallback near the boundaries.
     */
    private void p0Inter(double t, double[] result) {
        double idx = t / dt05;
        int lo = (int) Math.floor(idx);
        int last = 2 * Nt;
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
    private void solveL(int worker, double[] L, double startTime, double endTime) {
        if (endTime - startTime < EPS) return;

        double[] p0Buf = p0BufPool[worker];
        RungeKutta rk4 = rk4Pool[worker];

        RungeKutta.RhsFunction rhsInterp = (t, y, dydt) -> {
            p0Inter(t, p0Buf);
            LRhs(worker, y, p0Buf, dydt);
        };

        RungeKutta.RhsFunction rhsGrid = (t, y, dydt) -> {
            int idx = (int) Math.round(t / dt05);
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
                    rk4.step(t, boundary - t, L, L, NaTrunc + 1, rhsInterp);
                    clampNonNeg(L);
                    if (!rescaleL(worker, L)) { L[0] = Double.NaN; return; }
                    t = boundary;
                    tOnGrid = false;
                }
                epoch++;
                currentRates(worker, epoch);
            }

            if (tNext - t > EPS) {
                boolean nextOnGrid = (i <= lastIdx);
                boolean stepOnGrid = tOnGrid && nextOnGrid && Math.abs((tNext - t) - dt) < EPS;
                rk4.step(t, tNext - t, L, L, NaTrunc + 1, stepOnGrid ? rhsGrid : rhsInterp);
                clampNonNeg(L);
                if (!rescaleL(worker, L)) { L[0] = Double.NaN; return; }
                t = tNext;
                tOnGrid = nextOnGrid;
            }
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
     * Compute L over all external branches in one pass, grouping by tip height so the
     * cherry-sharing optimization is preserved within each sample-time group.
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
        double prevHeight = 0.0;

        for (int i = 0; i < numInvalid; i++) {
            double tipHeight = invalidTipHeights[i];
            double parentHeight = invalidTipParentHeights[i];
            int tipNum = invalidTipNums[i];

            if (tipHeight != groupTipHeight) {
                Arrays.fill(L, 0, NaTrunc + 1, 1.0);
                logScalePool[worker] = 0.0;
                groupTipHeight = tipHeight;
                prevHeight = tipHeight;
            }

            if (parentHeight != prevHeight) {
                solveL(worker, L, prevHeight, parentHeight);
                prevHeight = parentHeight;
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
        double endTime = tree.isRoot(root) ? originTime : tree.getNodeHeight(tree.getParent(root));

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

        if (originTime <= rootH) {
            if (DIAG_INF) System.err.println("ABD-Serial -Inf [rootHeight]: rootHeight="
                    + rootH + " >= originTime=" + originTime);
            return Double.NEGATIVE_INFINITY;
        }

        if (parametersDirty) {
            refreshRates();
            solveP0();
            rateStateDirty = true;
            parametersDirty = false;
        }

        modifiedNodeCount = 0;

        solveLTips();

        TreeUtils.postOrderTraversalList(tree, postOrder);

        if (taskPool == null) {
            for (int nodeNum : postOrder) {
                NodeRef node = tree.getNode(nodeNum);
                if (nodeValid[nodeNum] || tree.isExternal(node)
                        || (excludeRootBranch && tree.isRoot(node))) continue;
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
                    if (DIAG_INF) System.err.println("ABD-Serial -Inf [psi<=0]: psi=" + psi
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
        if (excludeRootBranch) {
            double rootHeight = tree.getNodeHeight(root);
            int left = tree.getChild(root, 0).getNumber();
            int right = tree.getChild(root, 1).getNumber();

            double leftL_0 = branchTopL[left][0];
            double rightL_0 = branchTopL[right][0];
            if (!Double.isFinite(leftL_0) || !Double.isFinite(rightL_0)
                    || leftL_0 <= 0.0 || rightL_0 <= 0.0) {
                if (DIAG_INF) System.err.println("ABD-Serial -Inf [rootChildL]: leftL_0="
                        + leftL_0 + " rightL_0=" + rightL_0 + " rootHeight=" + rootHeight);
                return Double.NEGATIVE_INFINITY;
            }

            double[] p0Buf = p0BufPool[0];
            p0Inter(rootHeight, p0Buf);

            double p0Root = p0Buf[0];
            if (!Double.isFinite(p0Root) || p0Root >= 1.0) {
                if (DIAG_INF) System.err.println("ABD-Serial -Inf [p0Root]: p0Root="
                        + p0Root + " rootHeight=" + rootHeight);
                return Double.NEGATIVE_INFINITY;
            }

            double logLPartial = Math.log(branchTopL[left][0]) + Math.log(branchTopL[right][0]);
            double logSurvival = -2.0 * Math.log(1.0 - p0Buf[0]);
            logLik += logLPartial + logSurvival + totalLogScale + logSamplingFactor;
        } else {
            double LRoot_0 = branchTopL[rootNum][0];
            if (!Double.isFinite(LRoot_0) || LRoot_0 <= 0.0) {
                if (DIAG_INF) System.err.println("ABD-Serial -Inf [rootL]: LRoot_0=" + LRoot_0);
                return Double.NEGATIVE_INFINITY;
            }

            double p0Origin = p0Grid[2 * Nt][0];
            if (!Double.isFinite(p0Origin) || p0Origin >= 1.0) {
                if (DIAG_INF) System.err.println("ABD-Serial -Inf [p0Origin]: p0Origin=" + p0Origin);
                return Double.NEGATIVE_INFINITY;
            }

            double logLPartial = Math.log(branchTopL[rootNum][0]);
            double logSurvival = -1.0 * Math.log(1.0 - p0Origin);
            logLik += logLPartial + logSurvival + totalLogScale + logSamplingFactor;
        }

        if (!Double.isFinite(logLik)) {
            if (DIAG_INF) System.err.println("ABD-Serial -Inf [nonFiniteLogLik]: logLik=" + logLik
                    + " totalLogScale=" + totalLogScale);
            return Double.NEGATIVE_INFINITY;
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
            if (excludeRootBranch && tree.isRoot(node)) continue;

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
     * legitimately depend on it:
     *   - the node itself (its branch was integrated against its parent height),
     *   - all its immediate children (their branchTopL is integrated up to this
     *     node's height; if the node moved or was rewired, those L curves
     *     belong to a different time interval now),
     *   - the full ancestor chain up to the root (each ancestor's branchTopL
     *     consumes the child L we just invalidated),
     *   - and each ancestor's other child along the chain, as a defense against
     *     topology operators that may not fire an event for every site.
     *
     * Wilson-Balding, SubtreeSlide, and exchange operators fire one
     * TreeChangedEvent per addChild call, so detach and reattach sites each
     * push their own event and the union of these per-event invalidations
     * covers the full set of dirty branches.
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
            int totalCols = 2 * Nt + 1;
            for (int i = 0; i < totalCols; i++) {
                System.arraycopy(p0Grid[i], 0, storedP0Grid[i], 0, NaTrunc + 1);
            }
            for (int k = 0; k < numEpochs; k++) {
                System.arraycopy(birthHaz[k], 0, storedBirthHaz[k], 0, Na + 1);
                System.arraycopy(deathHaz[k], 0, storedDeathHaz[k], 0, Na + 1);
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
                + ", p0Origin=" + p0Grid[2 * Nt][0] + ")\n";
    }

    public String toString() {
        return Double.toString(getLogLikelihood());
    }
}
