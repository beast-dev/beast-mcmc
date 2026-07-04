/*
 * BSplineBasisExpansion.java
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
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.math.SquaredCachedSplines;
import dr.math.matrix.SymmetricMatrixAccumulator;

/**
 * B-spline basis expansion wrapping {@link SquaredCachedSplines}.
 *
 * Canonical parameter layout:
 * <ul>
 *   <li>{@code includeIntercept = true}:  {@code [intercept, theta_0, ..., theta_{p-1}]}, dim = p+1</li>
 *   <li>{@code includeIntercept = false}: {@code [theta_0, ..., theta_{p-1}]}, dim = p</li>
 * </ul>
 *
 * This class does NOT call {@code addVariable(rateParameter)}: the likelihood owns
 * that registration (§6.1 of the refactor plan).  Sign canonicalization belongs
 * in the likelihood's {@code acceptState()}, not here (§14).
 *
 * @author Filippo Monti
 */
public final class BSplineBasisExpansion extends AbstractIntegratedBasisExpansion {

    private final SquaredCachedSplines splines;
    private final Parameter rateParameter;
    private final boolean includeIntercept;

    /**
     * @param splines          the backing spline model
     * @param rateParameter    the coefficient parameter in canonical layout (intercept first if used)
     * @param includeIntercept if true, psi(t) = [1, b_0(t), ..., b_{p-1}(t)]; otherwise psi(t) = [b_0..b_{p-1}]
     */
    public BSplineBasisExpansion(SquaredCachedSplines splines,
                                 Parameter rateParameter,
                                 boolean includeIntercept) {
        super("BSplineBasisExpansion");
        this.splines          = splines;
        this.rateParameter    = rateParameter;
        this.includeIntercept = includeIntercept;
        // Do NOT call addVariable(rateParameter): the likelihood owns that registration.
        // Add basis-shape sub-models here when adaptive knots are introduced.
    }

    @Override
    public int getDimension() {
        return rateParameter.getDimension();
    }

    @Override
    public Parameter getParameter() {
        return rateParameter;
    }

    @Override
    public void evaluateBasis(double t, BasisContext context, double[] buffer) {
        if (includeIntercept) {
            splines.evaluateAugmentedBasisInPlace(t, buffer); // [1, b_0, ..., b_{p-1}]
        } else {
            splines.evaluateBasisInPlace(t, buffer);          // [b_0, ..., b_{p-1}]
        }
    }

    /**
     * Returns eta(t) = gamma' psi(t) using the prebuilt polynomial from
     * {@link SquaredCachedSplines#evaluateSpline}, which is O(degree) after an amortized
     * O(spans*dim*degree) rebuild on the first call after a parameter change.  Much faster
     * than evaluating all basis functions individually as {@link #evaluateBasis} does.
     */
    @Override
    public double evaluateScalar(double t, BasisContext context) {
        double f = splines.evaluateSpline(t); // intercept + theta'b(t)
        if (includeIntercept) return f;
        return f - splines.getIntercept().getParameterValue(0); // theta'b(t) only
    }

    @Override
    public void addWeightedGramMatrix(double start, double end, double weight,
                                      BasisContext context, SymmetricMatrixAccumulator G) {
        if (includeIntercept) {
            splines.addScaledGramMatrix(start, end, weight, G);              // (p+1)×(p+1)
        } else {
            splines.addScaledGramMatrixNoIntercept(start, end, weight, G);   // p×p
        }
    }

    @Override
    public void integrateBasis(double start, double end, BasisContext context, double[] buffer) {
        if (includeIntercept) {
            splines.integrateAugmentedBasisInPlace(start, end, buffer); // [end-start, int b_0, ...]
        } else {
            splines.integrateBasisInPlace(start, end, buffer);          // [int b_0, ...]
        }
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index,
                                              Parameter.ChangeType type) {
        // Only basis-shape variables reach here; currently BSplineBasisExpansion has none.
        // When knot parameters are added they register here and propagate via fireModelChanged().
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
