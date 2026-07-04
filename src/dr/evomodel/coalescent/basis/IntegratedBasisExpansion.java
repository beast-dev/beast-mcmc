/*
 * IntegratedBasisExpansion.java
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
import dr.math.matrix.SymmetricMatrixAccumulator;

/**
 * A basis expansion {@code psi(t)} for the squared-basis coalescent.
 *
 * The rate function is {@code lambda(t) = epsilon + eta(t)^2}, where
 * {@code eta(t) = gamma' psi(t)} and {@code gamma} is the coefficient vector
 * returned by {@link #getParameter()}.
 *
 * The likelihood requires only two operations:
 * <ol>
 *   <li>{@link #evaluateBasis} — write {@code psi(t)} into a buffer;</li>
 *   <li>{@link #addWeightedGramMatrix} — accumulate {@code weight * integral psi psi' dt}.</li>
 * </ol>
 *
 * Implementations extend {@link dr.inference.model.AbstractModel}.
 * The likelihood wires listeners as follows (see §6.1 of the refactor plan):
 * <ul>
 *   <li>{@code addVariable(basis.getParameter())} — coefficient changes reach
 *       {@code handleVariableChangedEvent}; G is NOT invalidated.</li>
 *   <li>{@code addModel(basis)} — basis-shape changes fire
 *       {@code handleModelChangedEvent}; G IS invalidated.</li>
 * </ul>
 * Implementations must NOT call {@code addVariable(coefficients)} internally.
 *
 * @author Filippo Monti
 */
public interface IntegratedBasisExpansion extends Model {

    /** Dimension of {@code psi(t)} = dimension of {@link #getParameter()}. */
    int getDimension();

    /**
     * The BEAST {@link Parameter} for {@code gamma} in canonical layout
     * {@code [intercept, base_0, ..., base_{p-1}, covariates...]}.
     *
     * Registered directly by the likelihood via {@code addVariable(basis.getParameter())}.
     */
    Parameter getParameter();

    /**
     * Write {@code psi(t)} into {@code buffer[0..getDimension()-1]}.
     * No allocation.
     */
    void evaluateBasis(double t, BasisContext context, double[] buffer);

    /**
     * Returns {@code eta(t) = gamma' psi(t)}, the scalar whose square drives the rate.
     *
     * Implementations should use the most efficient available representation — e.g.
     * a prebuilt polynomial rather than evaluating all basis functions one by one.
     * {@link #evaluateBasis} is still needed for the gradient; this method serves
     * the likelihood evaluation path only.
     */
    double evaluateScalar(double t, BasisContext context);

    /**
     * Accumulate {@code weight * integral_start^end psi(t) psi(t)' dt} into {@code G}.
     * All indices are absolute; covariate wrappers write cross-terms at absolute offsets.
     */
    void addWeightedGramMatrix(double start, double end, double weight,
                               BasisContext context, SymmetricMatrixAccumulator G);

    /**
     * Write {@code integral_start^end psi(t) dt} into {@code buffer[0..getDimension()-1]}.
     * Used by covariate wrappers (ADDITIVE mode) to form cross-block contributions.
     * No allocation.
     */
    void integrateBasis(double start, double end, BasisContext context, double[] buffer);
}
