package dr.evomodel.speciation;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeUtils;
import dr.inference.model.*;
import dr.math.RungeKutta;
import dr.xml.Reportable;

/**
 * @author Frederik M. Andersen
 *
 * Compute likelihood for (a)symmetric time- and age-dependent birth-death model using PDE formulation
 * dp0/dt = dp0/da + mu(t,a) + lambda(t,a) * p0(t,a) * p0(t,0) - (lambda(t,a) + mu(t,a)) * p0(t,a)
 * dL/dt  = dL/da  + lambda(t,a) * (L(t,a)*p0(t,0) + p0(t,a)*L(t,0)) - (lambda(t,a) + mu(t,a)) * L(t,a)
 *
 * Partial likelihood PDE is solved over each branch in a post-order traversal by solving the system of
 * ODEs obtained from discretizing ages. ODEs are solved using a fixed-step Runge-Kutta method.
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
    private double da;
    private double dt;

    private final boolean symmetric;

    // Pre-allocated work arrays
    private final double[][] p0Grid;
    private final double[] p0t;
    private final double[][] branchTopL;

    private final int[] postOrder;
    private final double[] L0;
    private final double[] Lmerged;
    private double[] solveLCurr;
    private double[] solveLNext;

    private double logLikelihood = 0.0;
    private double storedLogLikelihood = 0.0;

    private final RungeKutta rk4;

    // Compute flags
    private boolean likelihoodKnown;
    private boolean storedLikelihoodKnown;

    public AgeDependentBirthDeathPDEModel(String name,
                                          Tree tree,
                                          Parameter birthScale,
                                          Parameter birthShape,
                                          Parameter deathScale,
                                          Parameter deathShape,
                                          Parameter epochTimes,
                                          int Na,
                                          int Nt,
                                          boolean symmetric) {
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
        addVariable(epochTimes);

        this.Na = Na;
        this.Nt = Nt;

        this.symmetric = symmetric;

        this.p0Grid = new double[Nt + 1][Na + 1];
        this.rk4 = new RungeKutta(Na + 1);
        this.p0t = new double[Na + 1];

        int totalNodes = tree.getNodeCount();
        this.branchTopL = new double[totalNodes][Na + 1];

        this.postOrder = new int[totalNodes];
        this.L0 = new double[Na + 1];
        this.Lmerged = new double[Na + 1];
        this.solveLCurr = new double[Na + 1];
        this.solveLNext = new double[Na + 1];
    }

    /**
     * Birth and death rate parametrization using linear-exponential age hazard.
     * rate(t, a) = scale(epoch) * (1 + b*a) * exp(-gamma*a)
     */
    private int getEpochIndex(double t) {
        for (int k = 0; k < epochTimes.getDimension() - 1; k++) {
            if (t < epochTimes.getParameterValue(k)) {
                return Math.min(k, birthScale.getDimension() - 1);
            }
        }
        return birthScale.getDimension() - 1;
    }

    protected double birthRate(double t, double a) {
        int epoch = getEpochIndex(t);
        double scale = birthScale.getParameterValue(epoch);
        double b = birthShape.getParameterValue(0);
        double gamma = birthShape.getParameterValue(1);
        return scale * (1.0 + b * a) * Math.exp(-gamma * a);
    }

    protected double deathRate(double t, double a) {
        int epoch = getEpochIndex(t);
        double scale = deathScale.getParameterValue(Math.min(epoch, deathScale.getDimension() - 1));
        double b = deathShape.getParameterValue(0);
        double gamma = deathShape.getParameterValue(1);
        return scale * (1.0 + b * a) * Math.exp(-gamma * a);
    }


    /**
     * 4th-order central finite difference for df/da (lower order near boundaries)
     */
    private double ageDeriv(double[] f, int j) {
        if (j >= 2 && j <= Na - 2) {
            // 4th-order central
            return (-f[j + 2] + 8.0 * f[j + 1] - 8.0 * f[j - 1] + f[j - 2]) / (12.0 * da);
        } else if (j == 1 || j == Na - 1) {
            // 2nd-order central
            return (f[j + 1] - f[j - 1]) / (2.0 * da);
        } else if (j == 0) {
            // 2nd-order forward
            return (-3.0 * f[0] + 4.0 * f[1] - f[2]) / (2.0 * da);
        } else { // j == Na
            // 2nd-order backward
            return (3.0 * f[Na] - 4.0 * f[Na - 1] + f[Na - 2]) / (2.0 * da);
        }
    }

    /**
     * Integrate p0 PDE from t=0 to t=originTime for all ages
     */
    private void p0Rhs(double t, double[] p0, double[] dp0dt) {
        double p0_0 = p0[0];

        for (int j = 0; j <= Na; j++) {
            double a = j * da;
            double lam = birthRate(t, a);
            double mu = deathRate(t, a);
            double r = lam + mu;

            dp0dt[j] = ageDeriv(p0, j) + mu + lam * p0[j] * p0_0 - r * p0[j];
        }
    }

    private void p0RhsSymmetric(double t, double[] p0, double[] dp0dt) {
        double p0_0 = p0[0];

        for (int j = 0; j <= Na; j++) {
            double a = j * da;
            double lam = birthRate(t, a);
            double mu = deathRate(t, a);
            double r = lam + mu;

            dp0dt[j] = ageDeriv(p0, j) + mu + lam * p0_0 * p0_0 - r * p0[j];
        }
    }

    private void solveP0() {
        // Boundary condition p0(0, a) = 0 for all a
        for (int j = 0; j <= Na; j++) {
            p0Grid[0][j] = 0.0;
        }

        double[] curr = new double[Na + 1];
        double[] next = new double[Na + 1];

        System.arraycopy(p0Grid[0], 0, curr, 0, Na + 1);

        RungeKutta.RhsFunction p0Rhs = symmetric ? this::p0RhsSymmetric : this::p0Rhs;

        for (int n = 0; n < Nt; n++) {
            double t = n * dt;
            rk4.step(t, dt, curr, next, Na + 1, p0Rhs);
            System.arraycopy(next, 0, p0Grid[n + 1], 0, Na + 1);

            // Swap
            double[] swap = curr;
            curr = next;
            next = swap;
        }
    }

    /**
     * Interpolate p0 at an arbitrary time
     */
    private void getP0AtTime(double t, double[] result) {
        double idx = t / dt;
        int lo = (int) Math.floor(idx);
        if (lo < 0) lo = 0;
        if (lo >= Nt) {
            System.arraycopy(p0Grid[Nt], 0, result, 0, Na + 1);
            return;
        }
        double frac = idx - lo;

        // Cubic interpolation using 4 points centered around lo..lo+1,
        // falling back to linear at the boundaries.
        if (lo >= 1 && lo + 2 <= Nt) {
            for (int j = 0; j <= Na; j++) {
                double fm1 = p0Grid[lo - 1][j];
                double f0  = p0Grid[lo][j];
                double f1  = p0Grid[lo + 1][j];
                double f2  = p0Grid[lo + 2][j];

                result[j] = f0 + 0.5 * frac * ((f1 - fm1) + frac * ((2.0 * fm1 - 5.0 * f0 + 4.0 * f1 - f2) + frac * ((-fm1 + 3.0 * f0 - 3.0 * f1 + f2))));
            }
        } else {
            for (int j = 0; j <= Na; j++) {
                result[j] = (1.0 - frac) * p0Grid[lo][j] + frac * p0Grid[lo + 1][j];
            }
        }
    }


    /**
     * Integrate L over a branch given boundary from offspring branches
     */
    private void lRhs(double t, double[] L, double[] dLdt) {
        getP0AtTime(t, p0t);
        double p0_0 = p0t[0];
        double L_0 = L[0];

        for (int j = 0; j <= Na; j++) {
            double a = j * da;
            double lam = birthRate(t, a);
            double r = lam + deathRate(t, a);

            dLdt[j] = ageDeriv(L, j) + lam * (L[j] * p0_0 + p0t[j] * L_0) - r * L[j];
        }
    }

    private void lRhsSymmetric(double t, double[] L, double[] dLdt) {
        getP0AtTime(t, p0t);
        double p0_0 = p0t[0];
        double L_0 = L[0];

        for (int j = 0; j <= Na; j++) {
            double a = j * da;
            double lam = birthRate(t, a);
            double r = lam + deathRate(t, a);

            dLdt[j] = ageDeriv(L, j) + 2.0 * lam * p0_0 * L_0 - r * L[j];
        }
    }

    private void solveL(double startTime, double endTime, double[] L0, double[] result) {
        double duration = endTime - startTime;
        if (duration <= 0) {
            System.arraycopy(L0, 0, result, 0, Na + 1);
            return;
        }

        int nSteps = Math.max(1, (int) Math.round(duration / dt));
        double localDt = duration / nSteps;

        System.arraycopy(L0, 0, solveLCurr, 0, Na + 1);

        RungeKutta.RhsFunction lRhs = symmetric ? this::lRhsSymmetric : this::lRhs;

        for (int n = 0; n < nSteps; n++) {
            double t = startTime + n * localDt;
            rk4.step(t, localDt, solveLCurr, solveLNext, Na + 1, lRhs);
            double[] swap = solveLCurr;
            solveLCurr = solveLNext;
            solveLNext = swap;
        }

        System.arraycopy(solveLCurr, 0, result, 0, Na + 1);
    }

    /**
     * Likelihood computation - Solving L PDE over each branch in a post-order direction
     */
    private double calculateLogLikelihood() {
        double t_or = epochTimes.getParameterValue(epochTimes.getDimension() - 1);
        this.da = t_or / Na;
        this.dt = t_or / Nt;

        solveP0();

        TreeUtils.postOrderTraversalList(tree, postOrder);

        for (int nodeNum : postOrder) {
            NodeRef node = tree.getNode(nodeNum);

            if (tree.isExternal(node)) {
                // Boundary L(0, a) = 1 for all a
                for (int j = 0; j <= Na; j++) L0[j] = 1.0;

                double nodeHeight = tree.getNodeHeight(node);
                double parentHeight;
                if (tree.isRoot(node)) {
                    parentHeight = epochTimes.getParameterValue(epochTimes.getDimension() - 1);
                } else {
                    parentHeight = tree.getNodeHeight(tree.getParent(node));
                }

                solveL(nodeHeight, parentHeight, L0, branchTopL[nodeNum]);

            } else if (tree.getChildCount(node) == 1) {
                // Degree-1 internal node (e.g. extra root in TreeModel): propagate L without merging
                NodeRef child = tree.getChild(node, 0);
                int childNum = child.getNumber();

                double nodeHeight = tree.getNodeHeight(node);
                double parentHeight;
                if (tree.isRoot(node)) {
                    parentHeight = epochTimes.getParameterValue(epochTimes.getDimension() - 1);
                } else {
                    parentHeight = tree.getNodeHeight(tree.getParent(node));
                }

                solveL(nodeHeight, parentHeight, branchTopL[childNum], branchTopL[nodeNum]);

            } else {
                // Boundary L_k(t_k, a) = lambda(t_k, a) * (L_i(t_k, a)*L_j(t_k, 0) + L_i(t_k, 0)*L_j(t_k, a))
                NodeRef left = tree.getChild(node, 0);
                NodeRef right = tree.getChild(node, 1);
                int leftNum = left.getNumber();
                int rightNum = right.getNumber();

                double Li_0 = branchTopL[leftNum][0];
                double Lj_0 = branchTopL[rightNum][0];

                double nodeHeight = tree.getNodeHeight(node);

                for (int j = 0; j <= Na; j++) {
                    double a = j * da;
                    double bRate = birthRate(nodeHeight, a);
                    if (symmetric) {
                        Lmerged[j] = bRate * Li_0 * Lj_0;
                    } else {
                        Lmerged[j] = bRate * (branchTopL[leftNum][j] * Lj_0 + Li_0 * branchTopL[rightNum][j]);
                    }
                }

                double parentHeight;
                if (tree.isRoot(node)) {
                    parentHeight = epochTimes.getParameterValue(epochTimes.getDimension() - 1);
                } else {
                    parentHeight = tree.getNodeHeight(tree.getParent(node));
                }

                solveL(nodeHeight, parentHeight, Lmerged, branchTopL[nodeNum]);
            }
        }

        NodeRef root = tree.getRoot();
        double likelihood = branchTopL[root.getNumber()][0];

        if (likelihood <= 0.0 || !Double.isFinite(likelihood)) {
            return Double.NEGATIVE_INFINITY;
        }

         int numSpeciations = tree.getInternalNodeCount();
        // Log and condition on non-extinction
        double p0Origin = p0Grid[Nt][0];
        if (p0Origin >= 1.0 || !Double.isFinite(p0Origin)) {
            return Double.NEGATIVE_INFINITY;
        }
        double logLik = Math.log(likelihood) - Math.log(1.0 - p0Origin);
        if (!symmetric) {
            logLik -= numSpeciations * Math.log(2.0);
        }
        if (!Double.isFinite(logLik)) {
            return Double.NEGATIVE_INFINITY;
        }
        return logLik;
    }

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
    }

    protected void handleModelChangedEvent(Model model, Object object, int index) {
        makeDirty();
    }

    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        makeDirty();
    }

    protected void storeState() {
        storedLogLikelihood = logLikelihood;
        storedLikelihoodKnown = likelihoodKnown;
    }

    protected void restoreState() {
        logLikelihood = storedLogLikelihood;
        likelihoodKnown = storedLikelihoodKnown;
    }

    protected void acceptState() {
    }

    public String getReport() {
        return "logLikelihood: " + getLogLikelihood();
    }

    public String toString() {
        return Double.toString(getLogLikelihood());
    }
}
