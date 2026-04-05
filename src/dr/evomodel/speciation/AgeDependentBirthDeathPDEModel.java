package dr.evomodel.speciation;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeUtils;
import dr.evomodel.tree.TreeChangedEvent;
import dr.inference.model.*;
import dr.math.RungeKutta;
import dr.xml.Reportable;
import java.util.Arrays;

/**
 * @author Frederik M. Andersen
 *
 * Compute likelihood for asymmetric time- and age-dependent birth-death model using PDE formulation
 *      dp0/dt = dp0/da + mu(t,a) + lambda(t,a) * p0(t,a) * p0(t,0) - (lambda(t,a) + mu(t,a)) * p0(t,a)
 *      dL/dt  = dL/da + lambda(t,a) * (L(t,a) * p0(t,0) + p0(t,a) * L(t,0)) - (lambda(t,a) + mu(t,a)) * L(t,a)
 *
 * or in the symmetric case
 *      dp0/dt = dp0/da + mu(t,a) + lambda(t,a) * p0(t,0)^2 - (lambda(t,a) + mu(t,a)) * p0(t,a)
 *      dL/dt  = dL/da + 2 * lambda(t,a) * L(t,0) * p0(t,0) - (lambda(t,a) + mu(t,a)) * L(t,a)
 *
 * Partial likelihood PDE is solved over each branch in a post-order traversal by solving the system of
 * ODEs obtained from discretizing ages. ODEs are solved using fixed-step RK4 with
 * h = min(da, dt) for CFL stability.
 *
 * Per-node caching: when only the tree changes (node height proposals), only affected nodes and
 * their ancestors are recomputed. The p0 grid (tree-independent) is reused from cache.
 */
public class AgeDependentBirthDeathPDEModel extends AbstractModelLikelihood implements Reportable {
    private final Tree tree;
    private final Parameter birthScale;
    private final Parameter birthShape;
    private final Parameter deathScale;
    private final Parameter deathShape;
    private final Parameter epochTimes;

    private final int Na;
    private final int Nt;
    private final double originTime;
    private final double da;
    private final double dt;
    private final double inv2da;
    private final double inv12da;
    private final double h;

    private final boolean symmetric;
    private final boolean useNodeCaching;

    // Pre-allocated work arrays
    private double[][] p0Grid;
    private double[][] storedP0Grid;
    private double[][] branchTopL;
    private double[][] storedBranchTopL;

    private final double[] p0t;
    private final int[] postOrder;
    private final double[] L0;
    private final double[] Lmerged;
    private double[] LCurr;
    private double[] LNext;
    private double[] p0Curr;
    private double[] p0Next;

    private double[] birthHaz;
    private double[] storedBirthHaz;
    private double[] deathHaz;
    private double[] storedDeathHaz;
    private double[] bScale;
    private double[] storedBScale;
    private double[] dScale;
    private double[] storedDScale;
    private double[] epochTimesC;
    private double[] storedCachedEpochTimes;
    private final int numEpochs;

    // Current epoch rates — set before each integration segment
    private double bScaleCurr;
    private double dScaleCurr;

    private double[] nodeLogScale;
    private double[] storedNodeLogScale;

    private double logLikelihood = 0.0;
    private double storedLogLikelihood = 0.0;

    // Compute flags
    private boolean likelihoodKnown = false;
    private boolean storedLikelihoodKnown = false;
    private boolean parametersDirty = true;
    private boolean storedParametersDirty = true;
    private boolean[] nodeValid;
    private boolean[] storedNodeValid;

    // Conditional store/restore tracking
    private int[] modifiedNodes;
    private int modifiedNodeCount = 0;
    private boolean rateStateDirty = false;
    // Set after restoreState uses pointer swaps — forces next storeState to full-copy
    private boolean storedStateDirty = true;

    private final RungeKutta rk45;

    public AgeDependentBirthDeathPDEModel(String name,
                                          Tree tree,
                                          Parameter birthScale,
                                          Parameter birthShape,
                                          Parameter deathScale,
                                          Parameter deathShape,
                                          Parameter epochTimes,
                                          double originTime,
                                          int Na,
                                          int Nt,
                                          boolean symmetric,
                                          boolean useNodeCaching) {
        super(name);

        this.tree = tree;
        this.birthScale = birthScale;
        this.birthShape = birthShape;
        this.deathScale = deathScale;
        this.deathShape = deathShape;
        this.epochTimes = epochTimes;

        if (tree instanceof Model) {
            addModel((Model) tree);
        }
        addVariable(birthScale);
        addVariable(birthShape);
        addVariable(deathScale);
        addVariable(deathShape);
        if (epochTimes != null) {
            addVariable(epochTimes);
        }

        this.Na = Na;
        this.Nt = Nt;

        this.originTime = originTime;
        this.da = originTime / Na;
        this.dt = originTime / Nt;
        this.inv2da = 1.0 / (2.0 * da);
        this.inv12da = 1.0 / (12.0 * da);
        this.h = Math.min(da, dt);

        this.symmetric = symmetric;
        this.useNodeCaching = useNodeCaching; // Should always be used, this is for testing

        this.p0Grid = new double[Nt + 1][Na + 1];
        this.storedP0Grid = new double[Nt + 1][Na + 1];
        this.rk45 = new RungeKutta(Na + 1);
        this.p0t = new double[Na + 1];

        int totalNodes = tree.getNodeCount();
        this.branchTopL = new double[totalNodes][Na + 1];
        this.storedBranchTopL = new double[totalNodes][Na + 1];

        this.postOrder = new int[totalNodes];
        this.L0 = new double[Na + 1];
        this.Lmerged = new double[Na + 1];
        this.LCurr = new double[Na + 1];
        this.LNext = new double[Na + 1];
        this.p0Curr = new double[Na + 1];
        this.p0Next = new double[Na + 1];

        this.nodeLogScale = new double[totalNodes];
        this.storedNodeLogScale = new double[totalNodes];
        this.nodeValid = new boolean[totalNodes];
        this.storedNodeValid = new boolean[totalNodes];
        this.modifiedNodes = new int[totalNodes];

        this.birthHaz = new double[Na + 1];
        this.storedBirthHaz = new double[Na + 1];
        this.deathHaz = new double[Na + 1];
        this.storedDeathHaz = new double[Na + 1];
        int numBoundaries = (epochTimes != null) ? epochTimes.getDimension() : 0;
        this.numEpochs = numBoundaries + 1;
        this.bScale = new double[numEpochs];
        this.storedBScale = new double[numEpochs];
        this.dScale = new double[numEpochs];
        this.storedDScale = new double[numEpochs];
        this.epochTimesC = new double[numBoundaries];
        this.storedCachedEpochTimes = new double[numBoundaries];
    }

    /**
     * Cache epoch scale values and age-hazard arrays
     * rate(t, a) = scale(t) * hazard(a)
     * with hazard(a) = (1 + b*a) * exp(-gamma*a)
     */
    private void cacheRates() {
        boolean constBirth = birthScale.getDimension() == 1;
        boolean constDeath = deathScale.getDimension() == 1;

        for (int k = 0; k < numEpochs; k++) {
            bScale[k] = birthScale.getParameterValue(constBirth ? 0 : k);
            dScale[k] = deathScale.getParameterValue(constDeath ? 0 : k);
        }
        for (int k = 0; k < epochTimesC.length; k++) {
            epochTimesC[k] = epochTimes.getParameterValue(k);
        }

        double bShapeLin = birthShape.getParameterValue(0);
        double bShapeExp = birthShape.getParameterValue(1);
        double dShapeLin = deathShape.getParameterValue(0);
        double dShapeExp = deathShape.getParameterValue(1);

        for (int j = 0; j <= Na; j++) {
            double a = j * da;
            birthHaz[j] = (1.0 + bShapeLin * a) * Math.exp(-bShapeExp * a);
            deathHaz[j] = (1.0 + dShapeLin * a) * Math.exp(-dShapeExp * a);
        }
    }

    //=================================
    // Solve p0 over [0, originTime]
    //=================================

    /*
     * Evaluate right hand side of p0 PDE
     */
    private void p0Rhs(double t, double[] p0, double[] dp0dt) {
        double p0_0 = p0[0];

        for (int j = 0; j <= Na; j++) {
            double deriv = ageDeriv(p0, j);
            double lam = bScaleCurr * birthHaz[j];
            double mu = dScaleCurr * deathHaz[j];
            double p0_b = symmetric ? p0_0 : p0[j];
            dp0dt[j] = deriv + mu + lam * p0_b * p0_0 - (lam + mu) * p0[j];
        }
    }

    /*
     * Fill p0Grid, i.e. solve p0(t, a) over the time and age grids
     * Using method of lines and a 4th order Runge-Kutta scheme
     */
    private void solveP0() {
        // Boundary condition p0(0, a) = 0
        for (int j = 0; j <= Na; j++) {
            p0Grid[0][j] = 0.0;
        }
        System.arraycopy(p0Grid[0], 0, p0Curr, 0, Na + 1);

        RungeKutta.RhsFunction p0Rhs = this::p0Rhs;

        double t = 0.0;
        int gridNext = 1;

        for (int epoch = 0; epoch < numEpochs; epoch++) {
            bScaleCurr = bScale[epoch];
            dScaleCurr = dScale[epoch];
            double epochEnd = (epoch < epochTimesC.length) ? epochTimesC[epoch] : originTime;

            while (t < epochEnd - 1e-15) {
                double gridEnd = gridNext * dt;
                double segEnd = Math.min(gridEnd, epochEnd);

                while (t < segEnd - 1e-15) {
                    double step = Math.min(h, segEnd - t);
                    rk45.step(t, step, p0Curr, p0Next, Na + 1, p0Rhs);
                    t += step;
                    double[] sw = p0Curr;
                    p0Curr = p0Next;
                    p0Next = sw;
                }

                // Store at grid point if we reached it
                if (t >= gridEnd - 1e-15) {
                    System.arraycopy(p0Curr, 0, p0Grid[gridNext], 0, Na + 1);
                    gridNext++;
                }
            }
        }
    }

    /**
     * Fourth order Catmull-Rom spline interpolation of p0, linear fallback at boundaries
     */
    private void p0Inter(double t, double[] result) {
        double idx = t / dt;
        int lo = (int) Math.floor(idx);
        if (lo < 0) lo = 0;
        if (lo >= Nt) {
            System.arraycopy(p0Grid[Nt], 0, result, 0, Na + 1);
            return;
        }
        double frac = idx - lo;

        if (lo >= 1 && lo + 2 <= Nt) {
            // Precompute Catmull-Rom basis weights
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

            for (int j = 0; j <= Na; j++) {
                result[j] = wm1 * rowM1[j] + w0 * row0[j] + w1 * row1[j] + w2 * row2[j];
            }
        } else {
            // Linear fallback near boundaries
            double w0 = 1.0 - frac;
            double[] row0 = p0Grid[lo];
            double[] row1 = p0Grid[lo + 1];
            for (int j = 0; j <= Na; j++) {
                result[j] = w0 * row0[j] + frac * row1[j];
            }
        }
    }

    /**
     * Interpolate p0 only with age 0
     */
    private double p0Inter0(double t) {
        double idx = t / dt;
        int lo = (int) Math.floor(idx);
        if (lo < 0) lo = 0;
        if (lo >= Nt) {
            return p0Grid[Nt][0];
        }
        double frac = idx - lo;

        if (lo >= 1 && lo + 2 <= Nt) {
            double frac2 = frac * frac;
            double frac3 = frac2 * frac;
            double wm1 = -0.5 * frac + frac2 - 0.5 * frac3;
            double w0  = 1.0 - 2.5 * frac2 + 1.5 * frac3;
            double w1  = 0.5 * frac + 2.0 * frac2 - 1.5 * frac3;
            double w2  = -0.5 * frac2 + 0.5 * frac3;
            return wm1 * p0Grid[lo - 1][0] + w0 * p0Grid[lo][0] + w1 * p0Grid[lo + 1][0] + w2 * p0Grid[lo + 2][0];
        } else {
            return (1.0 - frac) * p0Grid[lo][0] + frac * p0Grid[lo + 1][0];
        }
    }

    //====================================
    // Solve L over [startTime, endTime]
    //====================================

    /*
     * Evaluate right hand side of L PDE
     */
    private void lRhs(double t, double[] L, double[] dLdt) {
        double L_0 = L[0];
        double p0_0;
        if (symmetric) {
            p0_0 = p0Inter0(t);
        } else {
            p0Inter(t, p0t);
            p0_0 = p0t[0];
        }

        for (int j = 0; j <= Na; j++) {
            double deriv = ageDeriv(L, j);
            double lam = bScaleCurr * birthHaz[j];
            double r = lam + dScaleCurr * deathHaz[j];
            if (symmetric) {
                dLdt[j] = deriv + 2.0 * lam * p0_0 * L_0 - r * L[j];
            } else {
                dLdt[j] = deriv + lam * (L[j] * p0_0 + p0t[j] * L_0) - r * L[j];
            }
        }
    }

    /*
     * Solve L PDE using fourth order Runge-Kutta scheme
     */
    private void solveL(double startTime, double endTime, double[] boundaryL, double[] result) {
        double length = endTime - startTime;
        if (length <= 0) {
            System.arraycopy(boundaryL, 0, result, 0, Na + 1);
            return;
        }

        System.arraycopy(boundaryL, 0, LCurr, 0, Na + 1);

        RungeKutta.RhsFunction lRhs = this::lRhs;

        double t = startTime;
        int startEpoch = getEpochIndex(t);

        for (int epoch = startEpoch; epoch < numEpochs && t < endTime - 1e-15; epoch++) {
            bScaleCurr = bScale[epoch];
            dScaleCurr = dScale[epoch];

            double epochEnd = (epoch < epochTimesC.length) ? epochTimesC[epoch] : originTime;
            double segEnd = Math.min(endTime, epochEnd);

            while (t < segEnd - 1e-15) {
                double step = Math.min(h, segEnd - t);
                rk45.step(t, step, LCurr, LNext, Na + 1, lRhs);
                t += step;
                double[] swap = LCurr;
                LCurr = LNext;
                LNext = swap;
            }
        }

        System.arraycopy(LCurr, 0, result, 0, Na + 1);
    }

    /**
     * Rescale L array by L[0] to prevent underflow. Returns the local log-scale contribution.
     */
    private double rescaleL(double[] L) {
        double scale = L[0];
        if (scale > 0.0 && Double.isFinite(scale)) {
            double invScale = 1.0 / scale;
            for (int j = 0; j <= Na; j++) {
                L[j] *= invScale;
            }
            return Math.log(scale);
        }
        return 0.0;
    }


    // =========================
    // Likelihood computation
    // =========================

    private double calculateLogLikelihood() {
        if (originTime <= tree.getNodeHeight(tree.getRoot())) {
            return Double.NEGATIVE_INFINITY;
        }

        if (parametersDirty) {
            cacheRates();
            solveP0();
            rateStateDirty = true;
            parametersDirty = false;
        }

        modifiedNodeCount = 0;

        TreeUtils.postOrderTraversalList(tree, postOrder);

        for (int nodeNum : postOrder) {
            if (useNodeCaching && nodeValid[nodeNum]) {
                continue; // cached branchTopL and nodeLogScale are still valid
            }

            modifiedNodes[modifiedNodeCount] = nodeNum;
            modifiedNodeCount++;

            NodeRef node = tree.getNode(nodeNum);
            nodeLogScale[nodeNum] = 0.0;

            if (tree.isExternal(node)) {
                Arrays.fill(L0, 0, Na + 1, 1.0);

                double nodeHeight = tree.getNodeHeight(node);
                double parentHeight;
                if (tree.isRoot(node)) {
                    parentHeight = originTime;
                } else {
                    parentHeight = tree.getNodeHeight(tree.getParent(node));
                }

                solveL(nodeHeight, parentHeight, L0, branchTopL[nodeNum]);
                nodeLogScale[nodeNum] += rescaleL(branchTopL[nodeNum]);
            } else {
                NodeRef left = tree.getChild(node, 0);
                NodeRef right = tree.getChild(node, 1);
                int leftNum = left.getNumber();
                int rightNum = right.getNumber();

                double leftL_0 = branchTopL[leftNum][0];
                double rightL_0 = branchTopL[rightNum][0];

                double nodeHeight = tree.getNodeHeight(node);
                double bs = bScale[getEpochIndex(nodeHeight)];

                for (int j = 0; j <= Na; j++) {
                    double bRate = bs * birthHaz[j];
                    if (symmetric) {
                        Lmerged[j] = bRate * leftL_0 * rightL_0;
                    } else {
                        Lmerged[j] = bRate * (branchTopL[leftNum][j] * rightL_0 + leftL_0 * branchTopL[rightNum][j]);
                    }
                }

                nodeLogScale[nodeNum] += rescaleL(Lmerged);

                double parentHeight;
                if (tree.isRoot(node)) {
                    parentHeight = originTime;
                } else {
                    parentHeight = tree.getNodeHeight(tree.getParent(node));
                }

                solveL(nodeHeight, parentHeight, Lmerged, branchTopL[nodeNum]);
                nodeLogScale[nodeNum] += rescaleL(branchTopL[nodeNum]);
            }
        }

        Arrays.fill(nodeValid, true);

        double totalLogScale = 0.0;
        for (double v : nodeLogScale) {
            totalLogScale += v;
        }

        NodeRef root = tree.getRoot();
        double likelihood = branchTopL[root.getNumber()][0];

        if (likelihood <= 0.0 || !Double.isFinite(likelihood)) {
            return Double.NEGATIVE_INFINITY;
        }

        int numSpeciations = tree.getInternalNodeCount();
        double p0Origin = p0Grid[Nt][0];
        if (p0Origin >= 1.0 || !Double.isFinite(p0Origin)) {
            return Double.NEGATIVE_INFINITY;
        }
        double logLik = Math.log(likelihood) + totalLogScale - Math.log(1.0 - p0Origin);
        if (!symmetric) {
            logLik -= numSpeciations * Math.log(2.0);
        }

        if (!Double.isFinite(logLik)) {
            return Double.NEGATIVE_INFINITY;
        }
        return logLik;
    }



    // =======
    // Utils
    // =======

    /**
     * Find epoch index for a given time
     */
    private int getEpochIndex(double t) {
        int epoch = 0;
        while (epoch < epochTimesC.length && t >= epochTimesC[epoch]) {
            epoch++;
        }
        return epoch;
    }

    /**
     * Age finite difference (4th-order central, 2nd-order at boundaries) at a single index
     * Stencil: j=0 forward, j=1 central, j=2,...,Na-2 4th-order central, j=Na-1 central, j=Na backward.
     */
    private double ageDeriv(double[] f, int j) {
        if (j == 0) {
            return (-3.0 * f[0] + 4.0 * f[1] - f[2]) * inv2da;
        } else if (j == 1) {
            return (f[2] - f[0]) * inv2da;
        } else if (j <= Na - 2) {
            return (-f[j + 2] + 8.0 * f[j + 1] - 8.0 * f[j - 1] + f[j - 2]) * inv12da;
        } else if (j == Na - 1) {
            return (f[Na] - f[Na - 2]) * inv2da;
        } else {
            return (3.0 * f[Na] - 4.0 * f[Na - 1] + f[Na - 2]) * inv2da;
        }
    }

    /**
     * Invalidate a node, its children, and all ancestors up to the root
     */
    private void invalidateNode(NodeRef node) {
        for (int i = 0; i < tree.getChildCount(node); i++) {
            nodeValid[tree.getChild(node, i).getNumber()] = false;
        }

        NodeRef current = node;
        while (current != null) {
            nodeValid[current.getNumber()] = false;
            if (tree.isRoot(current)) break;
            current = tree.getParent(current);
        }
    }

    private void invalidateAllNodes() {
        Arrays.fill(nodeValid, false);
    }

    // ================
    // MCMC state management
    // ================

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

    public void makeDirty() {
        likelihoodKnown = false;
        parametersDirty = true;
        invalidateAllNodes();
    }

    protected void handleModelChangedEvent(Model model, Object object, int index) {
        if (model == tree) {
            if (useNodeCaching) {
                if (object instanceof TreeChangedEvent) {
                    TreeChangedEvent event = (TreeChangedEvent) object;
                    if (event.getNode() != null) {
                        invalidateNode(event.getNode());
                    } else {
                        invalidateAllNodes();
                    }
                }

                likelihoodKnown = false;
            } else {
                makeDirty();
            }
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

        if (storedStateDirty || rateStateDirty) {
            for (int i = 0; i <= Nt; i++) {
                System.arraycopy(p0Grid[i], 0, storedP0Grid[i], 0, Na + 1);
            }
            System.arraycopy(birthHaz, 0, storedBirthHaz, 0, birthHaz.length);
            System.arraycopy(deathHaz, 0, storedDeathHaz, 0, deathHaz.length);
            System.arraycopy(bScale, 0, storedBScale, 0, bScale.length);
            System.arraycopy(dScale, 0, storedDScale, 0, dScale.length);
            System.arraycopy(epochTimesC, 0, storedCachedEpochTimes, 0, epochTimesC.length);
            rateStateDirty = false;
        }

        if (storedStateDirty) {
            for (int i = 0; i < branchTopL.length; i++) {
                System.arraycopy(branchTopL[i], 0, storedBranchTopL[i], 0, Na + 1);
            }
            System.arraycopy(nodeLogScale, 0, storedNodeLogScale, 0, nodeLogScale.length);
            storedStateDirty = false;
        } else {
            for (int i = 0; i < modifiedNodeCount; i++) {
                int n = modifiedNodes[i];
                System.arraycopy(branchTopL[n], 0, storedBranchTopL[n], 0, Na + 1);
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

        if (rateStateDirty) {
            double[][] tmp2D;
            double[] tmpD;
            boolean[] tmpB;

            tmp2D = p0Grid; p0Grid = storedP0Grid; storedP0Grid = tmp2D;
            tmpD = birthHaz; birthHaz = storedBirthHaz; storedBirthHaz = tmpD;
            tmpD = deathHaz; deathHaz = storedDeathHaz; storedDeathHaz = tmpD;
            tmpD = bScale; bScale = storedBScale; storedBScale = tmpD;
            tmpD = dScale; dScale = storedDScale; storedDScale = tmpD;
            tmpD = epochTimesC; epochTimesC = storedCachedEpochTimes; storedCachedEpochTimes = tmpD;
            tmp2D = branchTopL; branchTopL = storedBranchTopL; storedBranchTopL = tmp2D;
            tmpD = nodeLogScale; nodeLogScale = storedNodeLogScale; storedNodeLogScale = tmpD;
            tmpB = nodeValid; nodeValid = storedNodeValid; storedNodeValid = tmpB;

            storedStateDirty = true;
        } else {
            for (int i = 0; i < modifiedNodeCount; i++) {
                int n = modifiedNodes[i];
                System.arraycopy(storedBranchTopL[n], 0, branchTopL[n], 0, Na + 1);
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
        return "logLikelihood: " + logLikelihood + "\n";
    }

    public String toString() {
        return Double.toString(getLogLikelihood());
    }
}
