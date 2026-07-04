/*
 * CovariateMode.java
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

/**
 * Controls how piecewise-constant covariates enter the basis expansion.
 *
 * @author Filippo Monti
 */
public enum CovariateMode {
    /**
     * Additive: {@code psi(t) = [base_psi(t), x(t)]}.
     * Dimension = {@code base.getDimension() + covariateDim}.
     * Coefficient vector: {@code [gamma(p), alpha(q)]}.
     */
    ADDITIVE,

    /**
     * Multiplicative (interaction): {@code psi(t) = [base_psi(t), x_1(t)*base_psi(t), ..., x_q(t)*base_psi(t)]}.
     * Dimension = {@code base.getDimension() * (1 + covariateDim)}.
     * Coefficient vector: {@code [gamma(p), delta_1(p), ..., delta_q(p)]}.
     * Each delta_k is the per-basis-function sensitivity to covariate k:
     *   eta(t) = (gamma + sum_k delta_k x_k(t))' psi(t).
     */
    MULTIPLICATIVE
}
