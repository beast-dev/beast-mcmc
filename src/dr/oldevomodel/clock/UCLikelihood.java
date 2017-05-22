/*
 * UCLikelihood.java
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

package dr.oldevomodel.clock;

import dr.evomodel.tree.TreeModel;
import dr.inference.model.Parameter;
import dr.math.MathUtils;
import dr.math.distributions.NormalDistribution;

/**
 * Calculates the likelihood of a set of rate changes in a tree, assuming that rates are lognormally distributed
 * cf Yang and Rannala 2006
 *
 * @author Michael Defoin Platel
 */
public class UCLikelihood extends RateEvolutionLikelihood {

    public UCLikelihood(TreeModel tree, Parameter ratesParameter, Parameter variance, Parameter rootRate, boolean isLogSpace) {

        super((isLogSpace) ? "LogNormally Distributed" : "Normally Distributed", tree, ratesParameter, rootRate, false);

        this.isLogSpace = isLogSpace;
        this.variance = variance;

        addVariable(variance);
    }

    /**
     * @return the log likelihood of the rate.
     */
    double branchRateChangeLogLikelihood(double foo1, double rate, double foo2) {
        double var = variance.getParameterValue(0);
        double meanRate = rootRateParameter.getParameterValue(0);


        if (isLogSpace) {
            final double logmeanRate = Math.log(meanRate);
            final double logRate = Math.log(rate);

            return NormalDistribution.logPdf(logRate, logmeanRate - (var / 2.), Math.sqrt(var)) - logRate;

        } else {
            return NormalDistribution.logPdf(rate, meanRate, Math.sqrt(var));
        }
    }

    double branchRateSample(double foo1, double foo2) {
        double meanRate = rootRateParameter.getParameterValue(0);

        double var = variance.getParameterValue(0);

        if (isLogSpace) {
            final double logMeanRate = Math.log(meanRate);

            return Math.exp(MathUtils.nextGaussian() * Math.sqrt(var) + logMeanRate - (var / 2.));
        } else {
            return MathUtils.nextGaussian() * Math.sqrt(var) + meanRate;
        }
    }

    private final Parameter variance;

    boolean isLogSpace = false;
}
