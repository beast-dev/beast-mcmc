/*
 * CompoundDiscontinuousPotentialProvider.java
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

import dr.inference.model.CompoundParameter;
import dr.inference.model.Parameter;

import java.util.ArrayList;
import java.util.List;

/**
 * Blockwise composition for discontinuous potential providers.
 *
 * @author Filippo Monti (powered by OpenAI)
 */
public class CompoundDiscontinuousPotentialProvider implements DiscontinuousPotentialProvider {

    private final List<DiscontinuousPotentialProvider> providers;
    private final CompoundParameter parameter;
    private final int[] providerByDimension;
    private final int[] localIndexByDimension;

    public CompoundDiscontinuousPotentialProvider(List<DiscontinuousPotentialProvider> providers) {
        if (providers == null || providers.isEmpty()) {
            throw new IllegalArgumentException("At least one discontinuous provider is required");
        }

        this.providers = new ArrayList<>(providers);
        this.parameter = new CompoundParameter("compoundDiscontinuous");

        int totalDimension = 0;
        for (DiscontinuousPotentialProvider provider : providers) {
            parameter.addParameter(provider.getParameter());
            totalDimension += provider.getDimension();
        }

        this.providerByDimension = new int[totalDimension];
        this.localIndexByDimension = new int[totalDimension];

        int offset = 0;
        for (int providerIndex = 0; providerIndex < providers.size(); providerIndex++) {
            DiscontinuousPotentialProvider provider = providers.get(providerIndex);
            for (int localIndex = 0; localIndex < provider.getDimension(); localIndex++) {
                providerByDimension[offset] = providerIndex;
                localIndexByDimension[offset] = localIndex;
                offset++;
            }
        }
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
    public void refresh() {
        for (DiscontinuousPotentialProvider provider : providers) {
            provider.refresh();
        }
    }

    @Override
    public void refreshAfterPositionUpdate() {
        for (DiscontinuousPotentialProvider provider : providers) {
            provider.refreshAfterPositionUpdate();
        }
    }

    @Override
    public double getLogDensity() {
        double total = 0.0;
        for (DiscontinuousPotentialProvider provider : providers) {
            double logDensity = provider.getLogDensity();
            if (Double.isInfinite(logDensity) && logDensity < 0.0) {
                return Double.NEGATIVE_INFINITY;
            }
            total += logDensity;
        }
        return total;
    }

    @Override
    public double getLogDensityAfterSingleCoordinateMove(int index, double proposedValue) {
        double total = 0.0;
        for (int providerIndex = 0; providerIndex < providers.size(); providerIndex++) {
            DiscontinuousPotentialProvider provider = providers.get(providerIndex);
            double contribution;
            if (providerByDimension[index] == providerIndex) {
                contribution = provider.getLogDensityAfterSingleCoordinateMove(localIndexByDimension[index], proposedValue);
            } else {
                contribution = provider.getLogDensity();
            }

            if (Double.isInfinite(contribution) && contribution < 0.0) {
                return Double.NEGATIVE_INFINITY;
            }
            total += contribution;
        }
        return total;
    }

    @Override
    public boolean isDiscontinuous(int index) {
        return providers.get(providerByDimension[index]).isDiscontinuous(localIndexByDimension[index]);
    }

    @Override
    public double getNextDiscontinuity(int index, double currentValue, double direction) {
        return providers.get(providerByDimension[index]).getNextDiscontinuity(
                localIndexByDimension[index], currentValue, direction);
    }

    @Override
    public double getPotentialDifferenceAcrossBoundary(int index, double boundary, double direction) {
        return providers.get(providerByDimension[index]).getPotentialDifferenceAcrossBoundary(
                localIndexByDimension[index], boundary, direction);
    }
}
