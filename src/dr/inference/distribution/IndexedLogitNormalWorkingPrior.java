package dr.inference.distribution;

import dr.inference.model.*;
import dr.math.MathUtils;
import dr.math.distributions.NormalDistribution;

import java.util.Arrays;

/**
 * Mixture prior for rho_b in (0,1).
 *
 * Let theta_b = logit(rho_b).
 *
 * If z_b = 0, then rho_b ~ Uniform(0,1).
 * If z_b = 1, then theta_b ~ Normal(logit(alpha_s), 1),
 * where s = atomAssignments[b].
 *
 * Therefore, for z_b = 1, the density in rho-space is obtained by
 * the change of variables:
 *
 *   p(rho_b) = p(theta_b) * |d theta_b / d rho_b|
 *
 * with |d theta_b / d rho_b| = 1 / (rho_b (1-rho_b)).
 */

/*
 * @author: Filippo Monti
*/

public class IndexedLogitNormalWorkingPrior extends AbstractModelLikelihood
        implements GradientProvider, ParametricMultivariateDistributionModel {

    public static final String ID = "indexedLogitNormalWorkingPrior";

    // Numerical safety
    private static final double EPS_RHO = 1e-12;   // clamp rho away from {0,1}
    private static final double EPS_A   = 1e-12;   // clamp alpha atoms away from {0,1}
    private static final double SIGMA = 1.0;

    // -----------------------
    // Model variables
    // -----------------------
    private final Parameter rho;                 // dimension n, in (0,1)
    private final Parameter alphaAtoms;          // dimension K, in (0,1)
    private final Parameter indicator;           // dimension n, 0/1-ish
    private final IndexedParameter atomAssignments;    // dimension n, deterministic index into alphaAtoms

    // -----------------------
    // Likelihood cache + state
    // -----------------------
    private boolean dirty = true;
    private double cachedLogL = Double.NaN;

    private boolean storedDirty = true;
    private double storedCachedLogL = Double.NaN;

    public IndexedLogitNormalWorkingPrior(final Parameter rho,
                                          final Parameter alphaAtoms,
                                          final Parameter indicator,
                                          final IndexedParameter atomAssignments) { // TODO generalize to allow transforms and distributions in the signature
        super(ID);

        this.rho = rho;
        this.alphaAtoms = alphaAtoms;
        this.indicator = indicator;
        this.atomAssignments = atomAssignments;
        final int n = rho.getDimension();

        if (indicator.getDimension() != n || atomAssignments.getDimension() != n) {
            throw new IllegalArgumentException("Dimension mismatch: parameter dim=" + n +
                    " indicator dim=" + indicator.getDimension() +
                    " atomAssignments dim=" + atomAssignments.getDimension());
        }
        addVariable(rho);
        addVariable(alphaAtoms);
        addVariable(indicator);

        addVariable(atomAssignments.getValuesParameter());
        addVariable(atomAssignments.getIndicesParameter());
    }

    // ============================================================
    // Likelihood (with caching)
    // ============================================================

    @Override
    public double getLogLikelihood() {
        if (dirty) {
            cachedLogL = computeLogLikelihood();
            dirty = false;
        }
        return cachedLogL;
    }

    private double computeLogLikelihood() {
        final int n = rho.getDimension();

        double sum = 0.0;

        for (int b = 0; b < n; b++) {
            final double rRaw = rho.getParameterValue(b);
            if (!Double.isFinite(rRaw)) return Double.NEGATIVE_INFINITY;

            if (!(rRaw >= 0.0 && rRaw <= 1.0)) {
                return Double.NEGATIVE_INFINITY;
            }
            final double r = clampOpen01(rRaw, EPS_RHO);

            if (isActive(b)) {
                sum += computeWorkingLikelihood(b, r);
            } else {
                sum += computeNeutralLikelihood(r);
            }
        }

        return sum;
    }

    private double computeNeutralLikelihood(double r) {
        return 0.0; // Uniform(0,1) log density
    }

    private double computeWorkingLikelihood(int b, double r) {
        final int s = toAtomIndex(atomAssignments.getParameterValue(b), alphaAtoms.getDimension());
        final double aRaw = alphaAtoms.getParameterValue(s);
        if (!Double.isFinite(aRaw)) return Double.NEGATIVE_INFINITY;
        final double a = clampOpen01(aRaw, EPS_A);
        double theta = logit(r);
        double mu = logit(a);
        double output = NormalDistribution.logPdf(theta, mu, SIGMA);
        output += -Math.log(r) - Math.log(1.0 - r); // Jacobian: log |d theta / d r| = -log r - log(1-r)
        return output;
    }

    private double gradientNeutralLikelihood(double r) {
        return 0.0; // Uniform(0,1) log density has zero gradient
    }

    private double gradientWorkingLikelihood(int b, double r) {
        final double theta = logit(r);
        final int s = toAtomIndex(atomAssignments.getParameterValue(b), alphaAtoms.getDimension());
        final double aRaw = alphaAtoms.getParameterValue(s);
        if (!Double.isFinite(aRaw)) {
            throw new IllegalArgumentException("Mean parameter must be in [0, 1]");
        }

        final double a = clampOpen01(aRaw, EPS_A);
        final double mu = logit(a);

        // d/dtheta log N(theta; mu, SIGMA) = -(theta - mu)/SIGMA^2
        final double dLogN_dTheta = NormalDistribution.gradLogPdf(theta, mu, SIGMA);

        // dtheta/dr = 1/(r(1-r))
        final double dTheta_dR = 1.0 / (r * (1.0 - r));

        // Jacobian term derivative: -1/r + 1/(1-r)
        final double dLogJac_dR = -1.0 / r + 1.0 / (1.0 - r);

        return dLogN_dTheta * dTheta_dR + dLogJac_dR;
    }


    // ============================================================
    // GradientProvider: gradient wrt rho
    // ============================================================

    @Override
    public int getDimension() {
        return rho.getDimension();
    }

    /**
     * Returns gradient w.r.t. rho coordinates.
     *
     * For each b:
     *  theta = logit(r)
     *  log p(r) = log N(theta; mu, 1) + log|dtheta/dr|
     *
     * d/dr log N(theta; mu, 1) = (d/dtheta log N) * (dtheta/dr)
     *   where d/dtheta log N = -(theta - mu)  (since SIGMA^2=1)
     *   and dtheta/dr = 1/(r(1-r))
     *
     * d/dr log|dtheta/dr| = d/dr[-log r - log(1-r)] = -1/r + 1/(1-r)
     */

    @Override
    public double[] getGradientLogDensity(final Object obj) {
        final int n = rho.getDimension();
        final double[] grad = new double[n];

        for (int b = 0; b < n; b++) {

            final double rRaw = rho.getParameterValue(b);
            if (!(rRaw >= 0.0 && rRaw <= 1.0) || !Double.isFinite(rRaw)) {
                throw new IllegalArgumentException("Parameter should be in [0, 1]");
            }

            final double r = clampOpen01(rRaw, EPS_RHO);

            if (isActive(b)) {
                grad[b] = gradientWorkingLikelihood(b, r);
            } else {
                grad[b] = gradientNeutralLikelihood(r);
            }
        }

        return grad;
    }

    // ============================================================
    // bookkeeping: dirtiness + MCMC state store/restore
    // ============================================================

    @Override
    public Model getModel() {
        return this;
    }

    @Override
    public void makeDirty() {
        dirty = true;
    }

    @Override
    public void handleModelChangedEvent(final Model model, final Object object, final int index) {
        makeDirty();
    }

    @Override
    public void handleVariableChangedEvent(final Variable variable,
                                           final int index,
                                           final Parameter.ChangeType type) {
        makeDirty();
    }

    private final boolean DEBUG_SR = false; // set to true to debug store/restore/accept calls

    @Override
    protected void storeState() {
        if (DEBUG_SR) System.err.println("STORE rhoWorkingPrior(LogitNormal)");
        storedDirty = dirty;
        storedCachedLogL = cachedLogL;
    }

    @Override
    protected void restoreState() {
        if (DEBUG_SR) System.err.println("RESTORE rhoWorkingPrior(LogitNormal)");
        dirty = storedDirty;
        cachedLogL = storedCachedLogL;

        if (DEBUG_SR) {
            final double recomputed = computeLogLikelihood();
            if (Double.isFinite(cachedLogL) && Double.isFinite(recomputed)) {
                if (Math.abs(cachedLogL - recomputed) > 1e-8) {
                    System.err.println("rhoWorkingPrior(LogitNormal) mismatch after restore: cached=" + cachedLogL + " recomputed=" + recomputed);
                    System.err.println("indicator=" + Arrays.toString(indicator.getParameterValues()));
                    System.err.println("rho=" + Arrays.toString(rho.getParameterValues()));
                    System.err.println("alphaAtoms=" + Arrays.toString(alphaAtoms.getParameterValues()));
                    System.err.println("atomAssignments.indices=" + Arrays.toString(atomAssignments.getIndicesParameter().getParameterValues()));
                }
            }
        }
    }

    @Override
    protected void acceptState() {
        storedDirty = dirty;
        storedCachedLogL = cachedLogL;
    }

    // ============================================================
    // ParametricMultivariateDistributionModel
    // ============================================================

    @Override
    public double[] nextRandom() {
        final int n = rho.getDimension();
        final double[] x = new double[n];

        if (indicator.getDimension() != n || atomAssignments.getDimension() != n) {
            Arrays.fill(x, 0.5);
            return x;
        }

        for (int b = 0; b < n; b++) {
            if (isActive(b)) {
                final double a = getAtomProbability(b);
                final double mu = logit(a);
                final double theta = NormalDistribution.quantile(MathUtils.nextDouble(), mu, SIGMA); // inverse-CDF sampling
                x[b] = logistic(theta);
            } else {
                x[b] = MathUtils.nextDouble();
            }
        }
        return x;
    }

    @Override
    public Variable<Double> getLocationVariable() {
        @SuppressWarnings("unchecked")
        final Variable<Double> v = (Variable<Double>) rho;
        return v;
    }

    /**
     * logPdf of an explicit vector x (in rho-space) under the current mixture assignments.
     */
    @Override
    public double logPdf(final double[] x) {
        final int n = rho.getDimension();
        if (x == null || x.length != n) return Double.NEGATIVE_INFINITY;
        if (indicator.getDimension() != n || atomAssignments.getDimension() != n) return Double.NEGATIVE_INFINITY;

        double sum = 0.0;
        for (int b = 0; b < n; b++) {
            final double rRaw = x[b];
            if (!(rRaw >= 0.0 && rRaw <= 1.0) || !Double.isFinite(rRaw)) return Double.NEGATIVE_INFINITY;

            final double r = clampOpen01(rRaw, EPS_RHO);
            final double theta = logit(r);

            if (isActive(b)) {
                final double a = getAtomProbability(b);
                final double mu = logit(a);
                sum += NormalDistribution.logPdf(theta, mu, SIGMA);
                sum += -Math.log(r) - Math.log(1.0 - r); // Jacobian
            }
        }

        return sum;
    }

    @Override
    public double[][] getScaleMatrix() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public double[] getMean() {
        final int n = rho.getDimension();
        final double[] mean = new double[n];

        for (int b = 0; b < n; b++) {
            if (isActive(b)) {
                final double a = getAtomProbability(b);
                mean[b] = logit(a);
            } else {
                mean[b] = 0.0;
            }
        }

        return mean;
    }
    @Override
    public String getType() {
        return ID;
    }

    // ============================================================
    // Helpers
    // ============================================================

    private boolean isActive(final int b) {
        return indicator.getParameterValue(b) >= 0.5;
    }

    private int getAtomIndex(final int b) {
        return toAtomIndex(atomAssignments.getParameterValue(b), alphaAtoms.getDimension());
    }

    private static int toAtomIndex(final double v, final int K) {
        final int s = (int) Math.rint(v);
        if (s < 0 || s >= K) {
            throw new IllegalArgumentException("atomAssignments out of range: " + s + " (K=" + K + ")");
        }
        return s;
    }

    private double getAtomProbability(final int b) {
        final double aRaw = alphaAtoms.getParameterValue(getAtomIndex(b));
        if (!Double.isFinite(aRaw)) {
            throw new IllegalArgumentException("Non-finite alpha atom at branch " + b);
        }
        return clampOpen01(aRaw, EPS_A);
    }

    private static double clampOpen01(final double v, final double eps) {
        if (v <= eps) return eps;
        final double hi = 1.0 - eps;
        if (v >= hi) return hi;
        return v;
    }

    private static double logit(final double p) {
        return Math.log(p / (1.0 - p));
    }

    private static double logistic(final double x) {
        if (x >= 0) {
            final double e = Math.exp(-x);
            return 1.0 / (1.0 + e);
        } else {
            final double e = Math.exp(x);
            return e / (1.0 + e);
        }
    }
}