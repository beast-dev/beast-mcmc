/*
 * DiscontinuousHmcUtils.java
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
 * Small shared helpers for the discontinuous HMC operators, factored out to
 * avoid copy-pasting them between {@link DiscontinuousHamiltonianMonteCarloOperator}
 * and {@link MixedDiscontinuousHamiltonianMonteCarloOperator}.
 *
 * @author Filippo Monti (powered by OpenAI)
 */
class DiscontinuousHmcUtils {

    private DiscontinuousHmcUtils() {
    }

    static void shuffleInPlace(int[] order) {
        for (int i = order.length - 1; i > 0; i--) {
            final int j = MathUtils.nextInt(i + 1);
            final int tmp = order[i];
            order[i] = order[j];
            order[j] = tmp;
        }
    }

    static void validatePositiveStepSize(double stepSize) {
        if (!(stepSize > 0.0)) {
            throw new IllegalArgumentException("Step size must be positive");
        }
    }

    static void validateRandomStepSizeFraction(double randomStepSizeFraction) {
        if (randomStepSizeFraction < 0.0 || randomStepSizeFraction > 1.0) {
            throw new IllegalArgumentException("Random step size fraction must be in [0, 1]");
        }
    }

    static void validatePositiveStepCount(int nSteps) {
        if (nSteps < 1) {
            throw new IllegalArgumentException("Number of steps must be at least one");
        }
    }

    static double drawStepSize(double stepSize, double randomStepSizeFraction) {
        if (randomStepSizeFraction == 0.0) {
            return stepSize;
        }
        return stepSize * (1.0 + randomStepSizeFraction * (MathUtils.nextDouble() - 0.5));
    }
}
