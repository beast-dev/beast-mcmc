/*
 * DiscontinuousMomentum.java
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

import dr.math.MathUtils;

/**
 * Laplace momentum helper for discontinuous coordinates.
 *
 * @author Filippo Monti (powered by OpenAI)
 */
public class DiscontinuousMomentum {

    private final double[] scales;

    public DiscontinuousMomentum(double[] scales) {
        if (scales == null || scales.length == 0) {
            throw new IllegalArgumentException("At least one scale is required");
        }

        this.scales = scales.clone();
        for (double scale : this.scales) {
            if (!(scale > 0.0)) {
                throw new IllegalArgumentException("All momentum scales must be positive");
            }
        }
    }

    public int getDimension() {
        return scales.length;
    }

    public double getScale(int index) {
        return scales[index];
    }

    public double[] drawMomentum() {
        final double[] momentum = new double[scales.length];
        for (int i = 0; i < momentum.length; i++) {
            momentum[i] = drawMomentum(i);
        }
        return momentum;
    }

    public double drawMomentum(int index) {
        final double sign = MathUtils.nextBoolean() ? 1.0 : -1.0;
        return sign * scales[index] * MathUtils.nextExponential(1.0);
    }

    public double getKineticEnergy(double momentum, int index) {
        return Math.abs(momentum) / scales[index];
    }

    public double getVelocity(double momentum, int index) {
        return sign(momentum) / scales[index];
    }

    public double updateMomentumAfterCrossing(double momentum, int index, double deltaU) {
        return momentum - scales[index] * sign(momentum) * deltaU;
    }

    public double reflectMomentum(double momentum) {
        return -momentum;
    }

    public double getTotalKineticEnergy(double[] momentum) {
        double total = 0.0;
        for (int i = 0; i < momentum.length; i++) {
            total += getKineticEnergy(momentum[i], i);
        }
        return total;
    }

    static double sign(double value) {
        if (value > 0.0) {
            return 1.0;
        } else if (value < 0.0) {
            return -1.0;
        }
        return 0.0;
    }
}
