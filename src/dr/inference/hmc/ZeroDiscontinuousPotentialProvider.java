/*
 * ZeroDiscontinuousPotentialProvider.java
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

import dr.inference.model.Parameter;

/**
 * Dimension-padding discontinuous potential with identically zero log density.
 *
 * Mixed discontinuous HMC uses one full parameter layout shared between smooth
 * gradients and discontinuous jump providers. Continuous-only blocks still need
 * a matching discontinuous provider block, but their local jump contribution is
 * zero because they are advanced by the smooth HMC flow.
 *
 * @author Filippo Monti (powered by OpenAI)
 */
public final class ZeroDiscontinuousPotentialProvider implements DiscontinuousPotentialProvider {

    private final Parameter parameter;

    public ZeroDiscontinuousPotentialProvider(final Parameter parameter) {
        if (parameter == null) {
            throw new IllegalArgumentException("parameter must be non-null");
        }
        this.parameter = parameter;
    }

    @Override
    public Parameter getParameter() {
        return parameter;
    }

    @Override
    public int getDimension() {
        return parameter.getDimension();
    }

    @Override
    public double getLogDensity() {
        return 0.0;
    }

    @Override
    public double getLogDensityAfterSingleCoordinateMove(final int index, final double proposedValue) {
        return 0.0;
    }

    @Override
    public boolean isDiscontinuous(final int index) {
        return false;
    }

    @Override
    public double getNextDiscontinuity(final int index,
                                       final double currentValue,
                                       final double direction) {
        return Double.NaN;
    }

    @Override
    public double getPotentialDifference(final int index,
                                         final double currentValue,
                                         final double proposedValue) {
        return 0.0;
    }
}
