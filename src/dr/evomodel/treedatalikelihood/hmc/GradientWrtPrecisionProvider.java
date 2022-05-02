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

import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.continuous.BranchSpecificGradient;
import dr.math.distributions.WishartSufficientStatistics;
import dr.math.interfaces.ConjugateWishartStatisticsProvider;

/**
 * @author Paul Bastide
 * @author Marc A. Suchard
 */

public interface GradientWrtPrecisionProvider {

//    double[] getGradientWrtPrecision(double[] vecV);

    double[] getGradientWrtPrecision(double[] vecV, double[] gradient);

    double[] getGradientWrtVariance(double[] vecP, double[] vecV, double[] gradient);

    ConjugateWishartStatisticsProvider getWishartStatistic();

    BranchSpecificGradient getBranchSpecificGradient();

    abstract class AbstractGradientWrtPrecisionProvider implements GradientWrtPrecisionProvider {
        int dim;

        public ConjugateWishartStatisticsProvider getWishartStatistic() {
            return null;
        }

        public BranchSpecificGradient getBranchSpecificGradient() {
            return null;
        }
    }

    class WishartGradientWrtPrecisionProvider extends AbstractGradientWrtPrecisionProvider {

        private final ConjugateWishartStatisticsProvider wishartStatistics;

        public WishartGradientWrtPrecisionProvider(ConjugateWishartStatisticsProvider wishartStatistics) {
            this.wishartStatistics = wishartStatistics;
            this.dim = wishartStatistics.getPrecisionParameter().getRowDimension();
        }

        public double[] getGradientWrtPrecision(double[] vecV, double[] gradient) {
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

        public double[] getGradientWrtVariance(double[] vecP, double[] vecV, double[] gradient) {
            MultivariateChainRule ruleI = new MultivariateChainRule.InverseGeneral(vecP);
            return ruleI.chainGradient(getGradientWrtPrecision(vecV, gradient));
        }

        @Override
        public ConjugateWishartStatisticsProvider getWishartStatistic() {
            return wishartStatistics;
        }
    }

    class BranchSpecificGradientWrtPrecisionProvider extends AbstractGradientWrtPrecisionProvider {

        private final BranchSpecificGradient branchSpecificGradient;

        public BranchSpecificGradientWrtPrecisionProvider(BranchSpecificGradient branchSpecificGradient) {
            this.branchSpecificGradient = branchSpecificGradient;
            this.dim = ((TreeDataLikelihood) branchSpecificGradient.getLikelihood()).getDataLikelihoodDelegate().getTraitDim();
        }

//        public double[] getGradientWrtPrecision(double[] vecV) {
//            double[] gradient = branchSpecificGradient.getGradientLogDensity(); // Get gradient wrt variance
//            return getGradientWrtPrecision(vecV, gradient);
//        }

        public double[] getGradientWrtPrecision(double[] vecV, double[] gradient) {
            MultivariateChainRule ruleI = new MultivariateChainRule.InverseGeneral(vecV);
            return ruleI.chainGradient(gradient);
        }

        public double[] getGradientWrtVariance(double[] vecP, double[] vecV, double[] gradient) {
            return gradient;
        }

        @Override
        public BranchSpecificGradient getBranchSpecificGradient() {
            return branchSpecificGradient;
        }
    }


}
