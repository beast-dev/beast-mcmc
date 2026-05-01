/*
 * CanonicalTipObservation.java
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

package dr.evomodel.treedatalikelihood.continuous.canonical;

import java.util.Arrays;

/**
 * Exact tip observation for the canonical OU tree pathway.
 *
 * <p>This type is intentionally not a {@link dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianState}.
 * Exact observations and missing coordinates are constraint-like factors rather than
 * ordinary finite normalized Gaussian states. Observed coordinates carry exact values;
 * unobserved coordinates are marked explicitly through {@link #observed}.
 */
public final class CanonicalTipObservation {

    public final int dim;
    public final double[] values;
    public final boolean[] observed;
    public int observedCount;

    public CanonicalTipObservation(final int dim) {
        if (dim < 1) {
            throw new IllegalArgumentException("dim must be positive");
        }
        this.dim = dim;
        this.values = new double[dim];
        this.observed = new boolean[dim];
        this.observedCount = 0;
    }

    public void setObserved(final double[] observedValues) {
        validateLength(observedValues, "observedValues");
        System.arraycopy(observedValues, 0, values, 0, dim);
        Arrays.fill(observed, true);
        observedCount = dim;
    }

    public void setMissing() {
        Arrays.fill(values, 0.0);
        Arrays.fill(observed, false);
        observedCount = 0;
    }

    public void setPartiallyObserved(final double[] sourceValues, final boolean[] observedMask) {
        validateLength(sourceValues, "sourceValues");
        if (observedMask.length != dim) {
            throw new IllegalArgumentException("observedMask dimension mismatch");
        }

        System.arraycopy(sourceValues, 0, values, 0, dim);
        int count = 0;
        for (int i = 0; i < dim; i++) {
            observed[i] = observedMask[i];
            if (observed[i]) {
                count++;
            }
        }
        observedCount = count;
    }

    public void copyFrom(final CanonicalTipObservation source) {
        if (source.dim != dim) {
            throw new IllegalArgumentException("Dimension mismatch: " + source.dim + " vs " + dim);
        }
        System.arraycopy(source.values, 0, values, 0, dim);
        System.arraycopy(source.observed, 0, observed, 0, dim);
        observedCount = source.observedCount;
    }

    private void validateLength(final double[] values, final String label) {
        if (values.length != dim) {
            throw new IllegalArgumentException(label + " dimension mismatch");
        }
    }
}
