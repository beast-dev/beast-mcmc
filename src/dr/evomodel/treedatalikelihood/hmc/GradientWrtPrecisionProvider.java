/*
 * GradientWrtPrecisionProvider.java
 *
 * Copyright (c) 2002-2018 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodel.treedatalikelihood.hmc;

import dr.evomodel.treedatalikelihood.continuous.BranchRateGradient;
import dr.math.distributions.WishartSufficientStatistics;
import dr.math.interfaces.ConjugateWishartStatisticsProvider;

/**
 * @author Paul Bastide
 * @author Marc A. Suchard
 */

public interface GradientWrtPrecisionProvider {

    double[] getGradientWrtPrecision(double[] vecV);

    ConjugateWishartStatisticsProvider getWishartStatistic();

    BranchRateGradient getBranchRateGradient();

    abstract class AbstractGradientWrtPrecisionProvider implements GradientWrtPrecisionProvider{
        int dim;

        public ConjugateWishartStatisticsProvider getWishartStatistic() {
            return null;
        }

        public BranchRateGradient getBranchRateGradient() {
            return null;
        }
    }

    class WishartGradientWrtPrecisionProvider extends AbstractGradientWrtPrecisionProvider {

        private final ConjugateWishartStatisticsProvider wishartStatistics;

        public WishartGradientWrtPrecisionProvider(ConjugateWishartStatisticsProvider wishartStatistics) {
            this.wishartStatistics = wishartStatistics;
            this.dim = wishartStatistics.getPrecisionParameter().getRowDimension();
        }

        public double[] getGradientWrtPrecision(double[] vecV) {
            // Statistics
            WishartSufficientStatistics wss = wishartStatistics.getWishartStatistics();
            double[] vecS = wss.getScaleMatrix();
            int numberTips = wss.getDf();

            return getGradientWrtPrecision(vecV, numberTips, vecS);
        }

        private double[] getGradientWrtPrecision(double[] vecV, int n, double[] vecS) {

            assert vecV.length == dim * dim;
            assert vecS.length == dim * dim;
            assert n > 0;

            double[] gradient = new double[dim * dim];

            for (int i = 0; i < dim * dim; ++i) {
                gradient[i] = 0.5 * (n * vecV[i] - vecS[i]);
            }

            return gradient;
        }

        @Override
        public ConjugateWishartStatisticsProvider getWishartStatistic() {
            return wishartStatistics;
        }
    }

    class BranchRateGradientWrtPrecisionProvider extends AbstractGradientWrtPrecisionProvider {

        private final BranchRateGradient branchRateGradient;

        public BranchRateGradientWrtPrecisionProvider(BranchRateGradient branchRateGradient) {
            this.branchRateGradient = branchRateGradient;
            this.dim = branchRateGradient.getDimension();
        }

        public double[] getGradientWrtPrecision(double[] vecV) {
            throw new RuntimeException("Not yet implemented");
        }

        @Override
        public BranchRateGradient getBranchRateGradient() {
            return branchRateGradient;
        }
    }


}
