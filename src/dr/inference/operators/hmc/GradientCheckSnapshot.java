/*
 * GradientCheckSnapshot.java
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

package dr.inference.operators.hmc;

import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.math.matrixAlgebra.ReadableVector;

/**
 * Captures a single gradient-check parameter snapshot and provides
 * deterministic restore + dirty behavior for analytic and numeric checks.
 */
final class GradientCheckSnapshot {

    private final Parameter parameter;
    private final Likelihood likelihood;
    private final double[] rawSnapshot;

    GradientCheckSnapshot(final Parameter parameter, final Likelihood likelihood) {
        this.parameter = parameter;
        this.likelihood = likelihood;
        this.rawSnapshot = parameter.getParameterValues();
    }

    double[] rawSnapshot() {
        return rawSnapshot.clone();
    }

    void setRawAndDirty(final double[] rawValues) {
        ReadableVector.Utils.setParameter(rawValues, parameter);
        likelihood.makeDirty();
    }

    void restoreAndDirty() {
        setRawAndDirty(rawSnapshot);
    }
}

