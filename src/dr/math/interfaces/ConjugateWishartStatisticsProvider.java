/*
 * ConjugateWishartStatisticsProvider.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.math.interfaces;

import dr.inference.model.MatrixParameterInterface;
import dr.math.distributions.WishartSufficientStatistics;

import java.util.List;

/**
 * An interface for classes that return conjugate outer products for Gibbs sampling of precision matrices
 *
 * @author Marc A. Suchard
 * @author Gabriela Cybis
 */

public interface ConjugateWishartStatisticsProvider {
    
    WishartSufficientStatistics getWishartStatistics();

    MatrixParameterInterface getPrecisionParameter();

    class CompoundWishartStatistics implements ConjugateWishartStatisticsProvider {

        private final List<ConjugateWishartStatisticsProvider> providers;
        private final int length;

        public CompoundWishartStatistics(List<ConjugateWishartStatisticsProvider> providers) {
            this.providers = providers;
            this.length = providers.get(0).getWishartStatistics().getScaleMatrix().length;
        }

        @Override
        public WishartSufficientStatistics getWishartStatistics() {

            int df = 0;
            double[] scale = new double[length];

            for (ConjugateWishartStatisticsProvider stat : providers) {
                df += stat.getWishartStatistics().getDf();
                double[] increment = stat.getWishartStatistics().getScaleMatrix();
                for (int i = 0; i < length; ++i) {
                    scale[i] += increment[i];
                }
            }

            return new WishartSufficientStatistics(df, scale);
        }

        @Override
        public MatrixParameterInterface getPrecisionParameter() {
            return providers.get(0).getPrecisionParameter();
        }
    };


}
