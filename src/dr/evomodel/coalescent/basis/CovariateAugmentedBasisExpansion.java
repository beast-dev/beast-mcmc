/*
 * CovariateAugmentedBasisExpansion.java
 *
 * Copyright (c) 2002-2026 Alexei Drummond, Andrew Rambaut and Marc Suchard
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
 *  BEAST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package dr.evomodel.coalescent.basis;

import dr.evomodel.coalescent.timeline.BasisContext;
import dr.inference.model.CompoundParameter;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.math.matrix.OffsetBlockSymmetricAccumulator;
import dr.math.matrix.SymmetricMatrixAccumulator;

/**
 * Decorator that augments any {@link IntegratedBasisExpansion} with piecewise-constant covariates.
 *
 * <h3>ADDITIVE mode</h3>
 * <pre>
 *   psi(t) = [base_psi(t), x(t)]     dim = p + q
 *   eta(t) = gamma' psi(t) + alpha' x(t)
 * </pre>
 * Gram:
 * <pre>
 *   [ integral h h'       (integral h) x' ]
 *   [ x (integral h)'     (b-a) x x'      ]
 * </pre>
 *
 * <h3>MULTIPLICATIVE mode (interaction)</h3>
 * <pre>
 *   psi(t) = [base_psi(t), x_1(t)*base_psi(t), ..., x_q(t)*base_psi(t)]   dim = p*(1+q)
 *   eta(t) = (gamma + sum_k delta_k x_k(t))' base_psi(t)
 * </pre>
 * Gram (since x is piecewise-constant so x is scalar within each segment):
 * <pre>
 *   block[0,   0  ] = G_base                   weight * integral psi psi'
 *   block[0,   k+1] = x_k * G_base             (cross)
 *   block[k+1, l+1] = x_k * x_l * G_base       (interaction, k<=l)
 * </pre>
 * The cross and interaction blocks are routed via {@link OffsetBlockSymmetricAccumulator}
 * (pre-allocated per-instance).
 *
 * Canonical parameter: {@code CompoundParameter([base.getParameter(), covariateCoefficients])}.
 * The likelihood registers the full compound via {@code addVariable(basis.getParameter())}.
 *
 * @author Filippo Monti
 */
public final class CovariateAugmentedBasisExpansion extends AbstractIntegratedBasisExpansion {

    private final IntegratedBasisExpansion base;
    private final Parameter covariateCoefficients;
    private final int covariateDim;
    private final CovariateMode mode;
    private final CompoundParameter parameter; // [base.getParameter(), covariateCoefficients]

    private final double[] baseIntegralBuffer; // workspace for base.integrateBasis
    private final double[] basePsiBuffer;      // workspace for evaluateBasis (MULTIPLICATIVE)

    public CovariateAugmentedBasisExpansion(IntegratedBasisExpansion base,
                                             Parameter covariateCoefficients,
                                             int covariateDim,
                                             CovariateMode mode) {
        super("CovariateAugmentedBasisExpansion");

        int p = base.getDimension();
        int expectedCoeffDim = (mode == CovariateMode.ADDITIVE) ? covariateDim : covariateDim * p;
        if (covariateCoefficients.getDimension() != expectedCoeffDim)
            throw new IllegalArgumentException("covariate coefficient dimension mismatch: " +
                    "expected " + expectedCoeffDim + ", got " + covariateCoefficients.getDimension());

        this.base                   = base;
        this.covariateCoefficients  = covariateCoefficients;
        this.covariateDim           = covariateDim;
        this.mode                   = mode;

        this.parameter = new CompoundParameter("augmentedBasisParameter");
        this.parameter.addParameter(base.getParameter());
        this.parameter.addParameter(covariateCoefficients);

        this.baseIntegralBuffer = new double[p];
        this.basePsiBuffer      = new double[p];

        // Register base as a sub-model so that basis-shape events propagate upward.
        // Do NOT call addVariable(covariateCoefficients): the likelihood owns that.
        addModel((Model) base);
    }

    @Override
    public int getDimension() {
        int p = base.getDimension();
        return (mode == CovariateMode.ADDITIVE) ? p + covariateDim : p * (1 + covariateDim);
    }

    @Override
    public Parameter getParameter() {
        return parameter;
    }

    @Override
    public void evaluateBasis(double t, BasisContext context, double[] buffer) {
        int p = base.getDimension();
        base.evaluateBasis(t, context, buffer); // fills buffer[0..p-1]
        double[] x = context.getCovariates();
        if (mode == CovariateMode.ADDITIVE) {
            for (int k = 0; k < covariateDim; k++) {
                buffer[p + k] = x[k];
            }
        } else { // MULTIPLICATIVE: buffer[p + k*p + i] = x[k] * psi_i(t)
            for (int k = 0; k < covariateDim; k++) {
                for (int i = 0; i < p; i++) {
                    buffer[p + k * p + i] = x[k] * buffer[i];
                }
            }
        }
    }

    @Override
    public double evaluateScalar(double t, BasisContext context) {
        double scalar = base.evaluateScalar(t, context); // gamma' psi(t)
        double[] x = context.getCovariates();
        if (mode == CovariateMode.ADDITIVE) {
            for (int k = 0; k < covariateDim; k++) {
                scalar += covariateCoefficients.getParameterValue(k) * x[k];
            }
        } else { // MULTIPLICATIVE: (gamma + sum_k delta_k x_k)' psi(t)
            int p = base.getDimension();
            base.evaluateBasis(t, context, basePsiBuffer);
            for (int k = 0; k < covariateDim; k++) {
                double xk = x[k];
                int offset = k * p; // delta_k starts at covariateCoefficients[offset]
                for (int i = 0; i < p; i++) {
                    scalar += covariateCoefficients.getParameterValue(offset + i) * xk * basePsiBuffer[i];
                }
            }
        }
        return scalar;
    }

    @Override
    public void addWeightedGramMatrix(double start, double end, double weight,
                                      BasisContext context, SymmetricMatrixAccumulator G) {
        int p = base.getDimension();
        double[] x = context.getCovariates();

        // Base block [0..p-1, 0..p-1]: always the same.
        base.addWeightedGramMatrix(start, end, weight, context, G);

        if (mode == CovariateMode.ADDITIVE) {
            // Cross block [0..p-1, p..p+q-1]: weight * (integral h) * x'.
            base.integrateBasis(start, end, context, baseIntegralBuffer);
            for (int i = 0; i < p; i++) {
                for (int k = 0; k < covariateDim; k++) {
                    G.add(i, p + k, weight * baseIntegralBuffer[i] * x[k]);
                }
            }
            // Bottom-right block [p..p+q-1, p..p+q-1]: weight * (end-start) * x x'.
            double length = end - start;
            for (int k = 0; k < covariateDim; k++) {
                for (int l = 0; l <= k; l++) {
                    G.add(p + k, p + l, weight * length * x[k] * x[l]);
                }
            }
        } else { // MULTIPLICATIVE
            // Cross block [0..p-1, p+k*p..p+(k+1)*p-1] = x_k * G_base.
            // Interaction block [p+k*p..p+(k+1)*p-1, p+l*p..p+(l+1)*p-1] = x_k*x_l * G_base.
            // Both routed via OffsetBlockSymmetricAccumulator so base.addWeightedGramMatrix
            // writes into the right block of G without copying.
            for (int k = 0; k < covariateDim; k++) {
                double xk = x[k];
                // Cross block (base × covariate k):
                base.addWeightedGramMatrix(start, end, weight * xk, context,
                        new OffsetBlockSymmetricAccumulator(G, 0, p + k * p));
                // Interaction blocks (covariate k × covariate l, k <= l):
                for (int l = k; l < covariateDim; l++) {
                    base.addWeightedGramMatrix(start, end, weight * xk * x[l], context,
                            new OffsetBlockSymmetricAccumulator(G, p + k * p, p + l * p));
                }
            }
        }
    }

    @Override
    public void integrateBasis(double start, double end, BasisContext context, double[] buffer) {
        int p = base.getDimension();
        base.integrateBasis(start, end, context, buffer); // fills buffer[0..p-1]
        double[] x = context.getCovariates();
        if (mode == CovariateMode.ADDITIVE) {
            double length = end - start;
            for (int k = 0; k < covariateDim; k++) {
                buffer[p + k] = x[k] * length;
            }
        } else { // MULTIPLICATIVE: buffer[p + k*p + i] = x[k] * integral psi_i dt
            for (int k = 0; k < covariateDim; k++) {
                double xk = x[k];
                for (int i = 0; i < p; i++) {
                    buffer[p + k * p + i] = xk * buffer[i];
                }
            }
        }
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index,
                                              Parameter.ChangeType type) {
        fireModelChanged();
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        fireModelChanged(object, index);
    }

    @Override protected void storeState()   {}
    @Override protected void restoreState() {}
    @Override protected void acceptState()  {}
}
