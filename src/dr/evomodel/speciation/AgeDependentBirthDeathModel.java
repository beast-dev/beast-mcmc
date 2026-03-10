package dr.evomodel.speciation;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.inference.model.*;

/**
 * @author Frederik M. Andersen
 *
 * Time- and Age-dependent birth-death model with spline-based rates on an NxN grid.
 *
 * Birth and death rates are:
 *   lambda(t, a) = exp(timeSpline(t) + ageSpline(a))
 *   mu(t, a)     = exp(timeSpline(t) + ageSpline(a))
 *
 * where timeSpline and ageSpline are B-spline expansions with additive level parameters.
 * Direct quadrature only (no FFT).
 */
public class AgeDependentBirthDeathModel extends AbstractModelLikelihood {

    private final Tree tree;

    // Spline parameters
    private final Parameter birthLevel;
    private final Parameter deathLevel;
    private final Parameter birthTimeCoefficients;
    private final Parameter birthAgeCoefficients;
    private final Parameter deathTimeCoefficients;
    private final Parameter deathAgeCoefficients;
    private final MatrixParameter timeBasis;
    private final MatrixParameter ageBasis;
    private final int numTimeBasis;
    private final int numAgeBasis;

    // Discretization
    private final int timeSteps;
    private final double h;

    // Cached spline evaluations
    private final double[] birthTimeSpline;
    private final double[] birthAgeSpline;
    private final double[] deathTimeSpline;
    private final double[] deathAgeSpline;
    private boolean splineCacheValid;

    // Cached rate grids
    private final double[][] birthRateGrid;
    private final double[][] deathRateGrid;
    private final double[][] expmR;

    // Extinction probability, branch survival probability and branch length densities
    private final double[] p0;
    private final double[] S;
    private final double[] branchDensBuffer;

    // State flags
    private boolean gridValid;
    private int sComputedUpTo;
    private boolean likelihoodKnown;
    private boolean storedLikelihoodKnown;
    private double logLikelihood;
    private double storedLogLikelihood;

    public AgeDependentBirthDeathModel(String name,
                                       Tree tree,
                                       Parameter birthLevel,
                                       Parameter deathLevel,
                                       Parameter birthTimeCoefficients,
                                       Parameter birthAgeCoefficients,
                                       Parameter deathTimeCoefficients,
                                       Parameter deathAgeCoefficients,
                                       MatrixParameter timeBasis,
                                       MatrixParameter ageBasis,
                                       double originTime,
                                       int timeSteps) {
        super(name);

        this.tree = tree;

        this.birthLevel = birthLevel;
        this.deathLevel = deathLevel;
        this.birthTimeCoefficients = birthTimeCoefficients;
        this.birthAgeCoefficients = birthAgeCoefficients;
        this.deathTimeCoefficients = deathTimeCoefficients;
        this.deathAgeCoefficients = deathAgeCoefficients;
        this.timeBasis = timeBasis;
        this.ageBasis = ageBasis;
        this.numTimeBasis = timeBasis.getColumnDimension();
        this.numAgeBasis = ageBasis.getColumnDimension();

        this.timeSteps = timeSteps;
        this.h = originTime / timeSteps;

        addVariable(birthLevel);
        addVariable(deathLevel);
        addVariable(birthTimeCoefficients);
        addVariable(birthAgeCoefficients);
        addVariable(deathTimeCoefficients);
        addVariable(deathAgeCoefficients);

        this.birthTimeSpline = new double[timeSteps + 1];
        this.birthAgeSpline = new double[timeSteps + 1];
        this.deathTimeSpline = new double[timeSteps + 1];
        this.deathAgeSpline = new double[timeSteps + 1];
        this.splineCacheValid = false;

        this.birthRateGrid = new double[timeSteps + 1][timeSteps + 1];
        this.deathRateGrid = new double[timeSteps + 1][timeSteps + 1];
        this.expmR = new double[timeSteps + 1][timeSteps + 1];

        this.p0 = new double[timeSteps + 1];
        this.S = new double[timeSteps + 1];
        this.branchDensBuffer = new double[timeSteps + 1];

        this.gridValid = false;
        this.likelihoodKnown = false;
    }

    private void evalSplines() {
        if (splineCacheValid) return;

        int nFreeTime = birthTimeCoefficients.getDimension();
        int nFreeAge = birthAgeCoefficients.getDimension();

        double bLevel = birthLevel.getParameterValue(0);
        double dLevel = deathLevel.getParameterValue(0);

        // Time splines (level absorbed here; last basis function dropped)
        for (int i = 0; i <= timeSteps; i++) {
            double birthSum = 0.0;
            double deathSum = 0.0;
            for (int k = 0; k < nFreeTime; k++) {
                birthSum += birthTimeCoefficients.getParameterValue(k) * timeBasis.getParameterValue(i, k);
                deathSum += deathTimeCoefficients.getParameterValue(k) * timeBasis.getParameterValue(i, k);
            }
            birthTimeSpline[i] = birthSum + bLevel;
            deathTimeSpline[i] = deathSum + dLevel;
        }

        // Age splines (last basis function dropped)
        for (int i = 0; i <= timeSteps; i++) {
            double birthSum = 0.0;
            double deathSum = 0.0;
            for (int l = 0; l < nFreeAge; l++) {
                birthSum += birthAgeCoefficients.getParameterValue(l) * ageBasis.getParameterValue(i, l);
                deathSum += deathAgeCoefficients.getParameterValue(l) * ageBasis.getParameterValue(i, l);
            }
            birthAgeSpline[i] = birthSum;
            deathAgeSpline[i] = deathSum;
        }

        splineCacheValid = true;
    }

    private double getBirthRate(int i, int j) {
        return Math.exp(birthTimeSpline[i] + birthAgeSpline[j]);
    }

    private double getDeathRate(int i, int j) {
        return Math.exp(deathTimeSpline[i] + deathAgeSpline[j]);
    }

    /*
     * Compute all rates, expmR, and extinction probabilities over the full grid.
     */
    private void computeGrid() {
        if (gridValid) return;

        evalSplines();

        // Base cases
        for (int i = 0; i <= timeSteps; i++) {
            birthRateGrid[i][0] = getBirthRate(i, 0);
            deathRateGrid[i][0] = getDeathRate(i, 0);
            expmR[i][0] = 1.0;
        }

        p0[0] = 0.0;
        S[0] = 1.0;
        sComputedUpTo = 0;

        // Compute rates, expmR, and p0 together for each age
        for (int m = 1; m <= timeSteps; m++) {
            // Rates and expmR at age m
            for (int i = 0; i <= timeSteps - m; i++) {
                birthRateGrid[i][m] = getBirthRate(i, m);
                deathRateGrid[i][m] = getDeathRate(i, m);
            }
            for (int i = timeSteps - m; i >= 0; i--) {
                double r1 = birthRateGrid[i][m] + deathRateGrid[i][m];
                double r2 = birthRateGrid[i + 1][m - 1] + deathRateGrid[i + 1][m - 1];
                expmR[i][m] = Math.exp(-h * (r1 + r2) / 2.0) * expmR[i + 1][m - 1];
            }

            // Extinction probability at index m
            double trap_sum_mu = 0.0;
            double trap_sum_lam = 0.0;

            trap_sum_mu += deathRateGrid[0][m] * expmR[0][m] / 2.0;

            for (int i = 1; i < m; i++) {
                int j = m - i;
                double expmR_ij = expmR[i][j];
                trap_sum_mu += deathRateGrid[i][j] * expmR_ij;
                trap_sum_lam += birthRateGrid[i][j] * expmR_ij * p0[i] * p0[i];
            }

            double expmR_m0 = expmR[m][0];
            trap_sum_mu += deathRateGrid[m][0] * expmR_m0 / 2.0;

            trap_sum_mu *= h;
            trap_sum_lam *= h;

            double a = h * birthRateGrid[m][0] * expmR_m0 / 2.0;
            double c = trap_sum_mu + trap_sum_lam;

            if (a == 0.0) {
                p0[m] = c;
            } else {
                double det = 1.0 - 4.0 * a * c;
                if (det < 0.0) {
                    p0[m] = 1.0 / (2.0 * a);
                } else {
                    p0[m] = (1.0 - Math.sqrt(det)) / (2.0 * a);
                }
            }
        }

        gridValid = true;
    }

    /*
     * Compute branch survival probability up to index m.
     */
    private void computeSUpTo(int idx) {
        if (idx <= sComputedUpTo) return;

        computeGrid();

        for (int i = sComputedUpTo + 1; i <= idx; i++) {
            double expmR_0idx = expmR[0][i];
            double trap_sum = 0.0;

            for (int j = 1; j < i; j++) {
                int k = i - j;
                trap_sum += birthRateGrid[j][k] * expmR[j][k] * p0[j] * S[j];
            }

            S[i] = (2.0 * h * trap_sum + expmR_0idx) /
                    (1.0 - h * birthRateGrid[i][0] * expmR[i][0] * p0[i]);
        }

        sComputedUpTo = idx;
    }

    /*
     * Compute branch length density for a node at height k with branch length l.
     */
    private double branchLengthDensity(int k, int l) {
        double p0_k = p0[k];
        double one_minus_p0_k = 1.0 - p0_k;

        branchDensBuffer[0] = birthRateGrid[k][0] * one_minus_p0_k * one_minus_p0_k;

        for (int m = 1; m <= l; m++) {
            double birthRate_k_m = birthRateGrid[k][m];
            double expmR_k_m = expmR[k][m];
            double trap_sum = birthRate_k_m * expmR_k_m * p0_k * branchDensBuffer[0] / 2.0;

            for (int r = 1; r < m; r++) {
                int i = k + r;
                int j = m - r;
                trap_sum += birthRateGrid[i][j] * expmR[i][j] * p0[i] * branchDensBuffer[r];
            }

            int i_end = k + m;
            double num = 2.0 * h * trap_sum + birthRate_k_m * expmR_k_m * one_minus_p0_k * one_minus_p0_k;
            double denom = 1.0 - h * birthRateGrid[i_end][0] * expmR[i_end][0] * p0[i_end];

            branchDensBuffer[m] = num / denom;
        }

        return Math.log(branchDensBuffer[l]) - Math.log(1.0 - p0[k + l]);
    }

    public double getLogLikelihood() {
        if (!likelihoodKnown) {
            logLikelihood = calculateLogLikelihood();
            likelihoodKnown = true;
        }
        return logLikelihood;
    }

    /*
     * Compute likelihood, simple loop over nodes
     */
    private double calculateLogLikelihood() {
        computeGrid();

        double logL = 0.0;

        for (int i = 0; i < tree.getInternalNodeCount(); i++) {
            NodeRef node = tree.getInternalNode(i);

            if (tree.isRoot(node)) {
                int k = (int) Math.round(tree.getNodeHeight(node) / h);
                logL += branchLengthDensity(k, timeSteps - k);
                continue;
            }

            int k = (int) Math.round(tree.getNodeHeight(node) / h);
            int l = (int) Math.round(tree.getBranchLength(node) / h);
            logL += branchLengthDensity(k, l);
        }

        for (int i = 0; i < tree.getExternalNodeCount(); i++) {
            NodeRef node = tree.getExternalNode(i);

            int l = (int) Math.round(tree.getBranchLength(node) / h);
            computeSUpTo(l);
            logL += Math.log(S[l]) - Math.log(1.0 - p0[l]);
        }

        if (Double.isNaN(logL) || logL == Double.POSITIVE_INFINITY) {
            return Double.NEGATIVE_INFINITY;
        }

        return logL;
    }

    protected void handleModelChangedEvent(Model model, Object object, int index) {
        likelihoodKnown = false;
        gridValid = false;
        splineCacheValid = false;
    }

    protected void storeState() {
        storedLikelihoodKnown = likelihoodKnown;
        storedLogLikelihood = logLikelihood;
    }

    protected void restoreState() {
        likelihoodKnown = storedLikelihoodKnown;
        logLikelihood = storedLogLikelihood;
        gridValid = false;
        splineCacheValid = false;
    }

    protected void acceptState() {
    }

    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        likelihoodKnown = false;
        gridValid = false;
        splineCacheValid = false;
    }

    public void makeDirty() {
        likelihoodKnown = false;
        gridValid = false;
        splineCacheValid = false;
    }

    public Model getModel() {
        return this;
    }

    public String toString() {
        return Double.toString(getLogLikelihood());
    }
}
