/*
 * CollapsedBridgeHamiltonianOperator.java
 *
 * Copyright © 2002-2026 the BEAST Development Team
 * http://beast.community/about
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * BEAST is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package dr.inference.operators.hmc;

import dr.inference.distribution.ExponentialTiltedStableDistribution;
import dr.inference.distribution.shrinkage.BayesianBridgeLikelihood;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.inference.operators.SimpleMCMCOperator;
import dr.math.MathUtils;

/**
 * Hamiltonian Monte Carlo on Bayesian-bridge coefficients theta with the local scales lambda integrated out,
 * followed by an exact redraw of every lambda from its exponentially tilted stable conditional.
 * <p>
 * The trajectory runs in z = log|theta| + kappa * sqrt|theta| with the signs of theta frozen, under the
 * collapsed density L(theta) * prod_k m(theta_k), m(theta) = exp(-|theta/tau|^alpha - theta^2 / (2 c^2)).
 * The target is the product form pi(theta, lambda) propto L * prod_k N(theta_k; 0, tau^2 lambda_k^2)
 * N(theta_k; 0, c^2) p0(lambda_k), the form the bridge Gibbs operator is exact for. The returned value
 * h = (K0 - K1) + (E1 - E0) - (P1 - P0) turns BEAST's acceptance into the collapsed-HMC ratio exp(-dH),
 * where E = sum_k [log m(a_k) + J(a_k)], J(a) = log a - log(1 + (kappa/2) sqrt a), and P is the score
 * terms (the prior) that BEAST adds itself. Fixed step size, no adaptation.
 * <p>
 * The score terms must be the joint minus the data likelihood. The move is valid only at chain temperature 1:
 * BEAST forms (newScore - oldScore) * temperature + h, so in a heated MC3 chain the subtracted P1 - P0 no longer
 * cancels. A tempered data likelihood with an untempered prior (power posterior) is fine.
 *
 * @author Marc A. Suchard
 */
public class CollapsedBridgeHamiltonianOperator extends SimpleMCMCOperator {

    public static final String OPERATOR_NAME = "collapsedBridgeHamiltonianOperator";

    public CollapsedBridgeHamiltonianOperator(Parameter rates,
                                              BayesianBridgeLikelihood bridge,
                                              GradientWrtParameterProvider likelihoodGradient,
                                              Likelihood scoreTerms,
                                              double[] massVariance,
                                              int nSteps, double stepSize, double kappa,
                                              double weight) {
        super(weight);
        this.rates = rates;
        this.bridge = bridge;
        this.localScale = bridge.getLocalScale();
        this.likelihoodGradient = likelihoodGradient;
        this.scoreTerms = scoreTerms;
        this.massVariance = massVariance;
        this.nSteps = nSteps;
        this.stepSize = stepSize;
        this.kappa = kappa;
        this.dim = rates.getDimension();
    }

    @Override
    public String getOperatorName() {
        return OPERATOR_NAME + "(" + rates.getParameterName() + ")";
    }

    @Override
    public double doOperation() {
        double[] momentum = new double[dim];
        for (int k = 0; k < dim; k++) {
            momentum[k] = MathUtils.nextGaussian() / Math.sqrt(massVariance[k]);
        }
        return proposeGivenMomentum(momentum);
    }

    /**
     * One collapsed-HMC proposal from the current state; leaves the final momentum in {@code momentum}.
     */
    double proposeGivenMomentum(double[] momentum) {

        readHyperparameters();

        final double[] theta = rates.getParameterValues();
        final double[] sign = new double[dim];
        final double[] z = new double[dim];
        for (int k = 0; k < dim; k++) {
            if (theta[k] == 0.0) {
                return Double.NEGATIVE_INFINITY;
            }
            sign[k] = Math.signum(theta[k]);
            z[k] = zOfA(Math.abs(theta[k]), kappa);
        }

        final double p0 = scoreTerms.getLogLikelihood();
        final double e0 = energy(theta);
        final double k0 = kineticEnergy(momentum);

        try {
            if (!leapFrog(z, momentum, sign)) {
                return Double.NEGATIVE_INFINITY;
            }
        } catch (ArithmeticException e) { // e.g. a failed eigendecomposition; rejected as in HamiltonianMonteCarloOperator
            return Double.NEGATIVE_INFINITY;
        }

        final double[] thetaNew = rates.getParameterValues();
        final double k1 = kineticEnergy(momentum);
        final double e1 = energy(thetaNew);

        for (int k = 0; k < dim; k++) {
            double omega = ExponentialTiltedStableDistribution.nextTiltedStable(
                    alpha / 2, Math.pow(thetaNew[k] / tau, 2));
            double lambda = Math.sqrt(1 / (2 * omega));
            if (!(lambda > 0.0 && lambda < Double.POSITIVE_INFINITY)) {
                return Double.NEGATIVE_INFINITY;
            }
            localScale.setParameterValueQuietly(k, lambda);
        }
        localScale.fireParameterChangedEvent();

        final double p1 = scoreTerms.getLogLikelihood();

        return (k0 - k1) + (e1 - e0) - (p1 - p0);
    }

    /**
     * nSteps leapfrog steps in (z, p); sets the rates at every position. Returns false on a non-finite value.
     */
    boolean leapFrog(double[] z, double[] momentum, double[] sign) {

        double[] gradient = gradZ(rates.getParameterValues());
        for (int k = 0; k < dim; k++) {
            momentum[k] += stepSize / 2 * gradient[k];
        }

        for (int step = 0; step < nSteps; step++) {

            for (int k = 0; k < dim; k++) {
                z[k] += stepSize * massVariance[k] * momentum[k];
                double a = aOfZ(z[k], kappa);
                if (!(a > 0.0 && a < Double.POSITIVE_INFINITY)) {
                    return false;
                }
                rates.setParameterValueQuietly(k, sign[k] * a);
            }
            rates.fireParameterChangedEvent();

            gradient = gradZ(rates.getParameterValues());
            final double scale = (step == nSteps - 1) ? stepSize / 2 : stepSize;
            for (int k = 0; k < dim; k++) {
                momentum[k] += scale * gradient[k];
                if (!(Math.abs(momentum[k]) < Double.POSITIVE_INFINITY)) {
                    return false;
                }
            }
        }
        return true;
    }

    /**
     * Gradient of log L + sum_k [log m(a_k) + J(a_k)] with respect to z; the rates must hold theta.
     */
    double[] gradZ(double[] theta) {

        final double[] logLikelihoodGradient = (likelihoodGradient != null) ?
                likelihoodGradient.getGradientLogDensity() : new double[dim];

        final double[] gradient = new double[dim];
        for (int k = 0; k < dim; k++) {
            final double a = Math.abs(theta[k]);
            final double r = Math.sqrt(a);
            final double fp = 1 + kappa / 2 * r;
            gradient[k] = (theta[k] * logLikelihoodGradient[k] + 1 - alpha * Math.pow(a / tau, alpha)
                    - theta[k] * theta[k] * inverseSlabSquared - (kappa / 4) * r / fp) / fp;
        }
        return gradient;
    }

    /**
     * sum_k [log m(a_k) + J(a_k)].
     */
    double energy(double[] theta) {
        double sum = 0.0;
        for (int k = 0; k < dim; k++) {
            final double a = Math.abs(theta[k]);
            sum += -Math.pow(a / tau, alpha) - a * a * inverseSlabSquared / 2
                    + Math.log(a) - Math.log1p(kappa / 2 * Math.sqrt(a));
        }
        return sum;
    }

    private double kineticEnergy(double[] momentum) {
        double sum = 0.0;
        for (int k = 0; k < dim; k++) {
            sum += massVariance[k] * momentum[k] * momentum[k];
        }
        return sum / 2;
    }

    void readHyperparameters() {
        tau = bridge.getGlobalScale().getParameterValue(0);
        alpha = bridge.getExponent().getParameterValue(0);
        Parameter slab = bridge.getSlabWidth();
        inverseSlabSquared = (slab != null) ? 1 / Math.pow(slab.getParameterValue(0), 2) : 0.0;
    }

    static double zOfA(double a, double kappa) {
        return Math.log(a) + kappa * Math.sqrt(a);
    }

    /**
     * Inverse of z = log a + kappa sqrt(a): a = (2/kappa)^2 W0((kappa/2) e^(z/2))^2. Newton on
     * g(u) = 2u + kappa e^u - z (u = log sqrt a; convex and increasing, started to the right of the root, so it
     * converges monotonically), then Newton on f(s) = 2 log s + kappa s - z until the step is at most 1 ulp.
     * Returns NaN if either stage does not converge.
     */
    static double aOfZ(double z, double kappa) {

        if (kappa == 0.0) {
            return Math.exp(z);
        }
        if (!(Math.abs(z) < Double.POSITIVE_INFINITY)) {
            return Double.NaN;
        }

        double u = z / 2; // g(z / 2) = kappa e^(z/2) > 0
        if (z > kappa) {
            u = Math.min(u, Math.log(z / kappa)); // g(log(z / kappa)) = 2 log(z / kappa) > 0
        }
        boolean converged = false;
        for (int i = 0; i < MAX_NEWTON_ITERATIONS && !converged; i++) {
            final double e = kappa * Math.exp(u);
            final double du = (2 * u + e - z) / (2 + e);
            u -= du;
            converged = Math.abs(du) <= 1E-8 * Math.max(1.0, Math.abs(u));
        }
        if (!converged) {
            return Double.NaN;
        }

        double s = Math.exp(u);
        for (int i = 0; i < MAX_NEWTON_ITERATIONS; i++) {
            final double f = newtonFunction(s, kappa, z);
            final double ds = f * s / (2 + kappa * s);
            if (Math.abs(ds) <= Math.ulp(s)) {
                final double square = s * s; // (s - ds)^2 = s^2 - 2 s ds, with s^2 exact and one final rounding
                return square + (productError(s, s, square) - square * (2 * f / (2 + kappa * s)));
            }
            s -= ds;
        }
        return Double.NaN;
    }

    /**
     * f(s) = 2 log s + kappa s - z, accurate enough for Newton to reach a step of 1 ulp of s: log s = e ln2 + log m
     * (s = m 2^e, ln2 split into an exactly multipliable high part and a low part) and kappa s = product + its exact
     * rounding error, so that the large terms cancel first. Plain 2 log s + kappa s - z is off by several ulps.
     */
    private static double newtonFunction(double s, double kappa, double z) {
        int e = Math.getExponent(s);
        double m = Math.scalb(s, -e);
        if (m > SQRT_TWO) {
            m /= 2;
            e += 1;
        }
        final double product = kappa * s;
        return ((2 * e * LN2_HI - z) + product) + (productError(kappa, s, product) + 2 * (e * LN2_LO + Math.log(m)));
    }

    private static double productError(double x, double y, double product) { // x * y - product, exactly (Dekker)
        final double xs = SPLIT * x, xh = xs - (xs - x), xl = x - xh;
        final double ys = SPLIT * y, yh = ys - (ys - y), yl = y - yh;
        return ((xh * yh - product) + xh * yl + xl * yh) + xl * yl;
    }

    private static final int MAX_NEWTON_ITERATIONS = 100;
    private static final double LN2_HI = 0x1.62e42feep-1;          // fdlibm ln2_hi: 32 significant bits
    private static final double LN2_LO = 0x1.a39ef35793c76p-33;    // ln 2 - LN2_HI
    private static final double SQRT_TWO = Math.sqrt(2.0);
    private static final double SPLIT = 134217729.0;               // 2^27 + 1

    private final Parameter rates;
    private final BayesianBridgeLikelihood bridge;
    private final Parameter localScale;
    private final GradientWrtParameterProvider likelihoodGradient;
    private final Likelihood scoreTerms;
    private final double[] massVariance;
    private final int nSteps;
    private final double stepSize;
    private final double kappa;
    private final int dim;

    private double tau;
    private double alpha;
    private double inverseSlabSquared;
}
