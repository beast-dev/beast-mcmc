/*
 * DiagonalPriorPreconditioner.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
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
 *  BEAST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 *
 */

package dr.inference.hmc;

import dr.inference.model.PriorPreconditioningProvider;
import dr.math.matrixAlgebra.Vector;
import dr.xml.Reportable;

import java.util.Arrays;

/**
 * Fixed vector of prior standard deviations for diagonal HMC preconditioning.
 *
 * @author Filippo Monti
 */
public class DiagonalPriorPreconditioner implements PriorPreconditioningProvider, Reportable {

    private final String id;
    private final double[] standardDeviations;

    public DiagonalPriorPreconditioner(final String id, final double[] standardDeviations) {
        if (standardDeviations.length == 0) {
            throw new IllegalArgumentException("Prior preconditioner requires at least one standard deviation");
        }
        this.id = id;
        this.standardDeviations = Arrays.copyOf(standardDeviations, standardDeviations.length);
        for (int i = 0; i < this.standardDeviations.length; i++) {
            final double sd = this.standardDeviations[i];
            if (!Double.isFinite(sd) || sd <= 0.0) {
                throw new IllegalArgumentException(
                        "Prior preconditioner standard deviation must be positive and finite at index " + i);
            }
        }
    }

    @Override
    public double getStandardDeviation(final int index) {
        return standardDeviations[index];
    }

    @Override
    public int getDimension() {
        return standardDeviations.length;
    }

    @Override
    public String getReport() {
        final StringBuilder sb = new StringBuilder("diagonalPriorPreconditioner Report\n\n");
        sb.append("id: ").append(id).append("\n\n");
        sb.append("dimension: ").append(getDimension()).append("\n\n");
        sb.append("Prior SDs: ");
        sb.append(new Vector(standardDeviations));
        sb.append("\n\n");
        return sb.toString();
    }
}
