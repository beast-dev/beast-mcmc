/*
 * TruncatedSoftBoundCauchyDistribution.java
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

package dr.inference.distribution;

import dr.inference.model.GradientProvider;
import dr.math.UnivariateFunction;
import dr.math.distributions.Distribution;
import dr.util.Author;
import dr.util.Citable;
import dr.util.Citation;
import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.AbstractContinuousDistribution;

import java.util.Collections;
import java.util.List;

/**
 * The uniform distribution with soft bounds as Eq. 17 in Yang and Rannala, 2006
 * reference: paml4.8 package by Dr. Ziheng Yang
 *
 * @author Xiang Ji
 * @author Marc Suchard
 */
public class SoftBoundUniformDistribution extends AbstractContinuousDistribution implements Distribution, GradientProvider, Citable {

    private final double lowerBound;
    private final double upperBound;
    private final double lowerProbability;
    private final double upperProbability;

    public SoftBoundUniformDistribution (double lowerBound,
                                         double upperBound,
                                         double lowerProbability,
                                         double upperProbability) {

        this.lowerBound = lowerBound;
        this.upperBound = upperBound;
        this.lowerProbability = lowerProbability;
        this.upperProbability = upperProbability;
    }

    @Override
    public double pdf(double x) {
        if (x > lowerBound && x < upperBound) {
            return (1.0 - lowerProbability - upperProbability) / (upperBound - lowerBound);
        } else if (x < lowerBound) {
            final double thetaL = (1.0 - lowerProbability - upperProbability)*lowerBound/(lowerProbability *(upperBound-lowerBound));
            return lowerProbability *thetaL/lowerBound * Math.pow(x/lowerBound, thetaL-1);
        } else {
            final double thetaR = (1.0 - lowerProbability - upperProbability) / (upperProbability * (upperBound - lowerBound));
            return upperProbability * thetaR / Math.exp(thetaR * (x - upperBound));
        }
    }

    /*
    Reference: line 2290 - line 2303 of mcmctree.c from paml4.8 package by Dr. Ziheng Yang
     */
    @Override
    public double logPdf(double x) {
        if (x > lowerBound && x < upperBound) {
            return Math.log((1.0 - lowerProbability - upperProbability) / (upperBound - lowerBound));
        } else if (x < lowerBound) {
            final double thetaL = (1.0 - lowerProbability - upperProbability)*lowerBound/(lowerProbability *(upperBound-lowerBound));
            return Math.log(lowerProbability *thetaL/lowerBound) + (thetaL-1)*Math.log(x/lowerBound);
        } else {
            final double thetaR = (1.0 - lowerProbability - upperProbability) / (upperProbability * (upperBound - lowerBound));
            return Math.log(upperProbability * thetaR) - thetaR * (x - upperBound);
        }
    }

    @Override
    public double cdf(double x) {
        throw new RuntimeException("Not yet implemented!");
    }

    @Override
    public double quantile(double y) {
        throw new RuntimeException("Not yet implemented!");
    }

    @Override
    public double mean() {
        throw new RuntimeException("Not yet implemented!");
    }

    @Override
    public double variance() {
        throw new RuntimeException("Not yet implemented!");
    }

    @Override
    public UnivariateFunction getProbabilityDensityFunction() {
        throw new RuntimeException("Not yet implemented!");
    }

    @Override
    protected double getInitialDomain(double v) {
        throw new RuntimeException("Not yet implemented!");
    }

    @Override
    protected double getDomainLowerBound(double v) {
        throw new RuntimeException("Not yet implemented!");
    }

    @Override
    protected double getDomainUpperBound(double v) {
        throw new RuntimeException("Not yet implemented!");
    }

    @Override
    public double cumulativeProbability(double v) throws MathException {
        throw new RuntimeException("Not yet implemented!");
    }

    @Override
    public Citation.Category getCategory() {
        return Citation.Category.PRIOR_MODELS;
    }

    @Override
    public String getDescription() {
        return "Uniform distribution with soft bound.";
    }

    @Override
    public List<Citation> getCitations() {
        return Collections.singletonList(CITATION);
    }

    public static Citation CITATION = new Citation(
            new Author[]{
                    new Author("Z", "Yang"),
                    new Author("B", "Rannala")
            },
            "Bayesian Estimation of Species Divergence Times Under a Molecular Clock Using Multiple Fossil Calibrations with Soft Bounds",
            2006,
            "Molecular Biology and Evolution",
            23,
            212,
            226,
            Citation.Status.PUBLISHED
    );

    @Override
    public int getDimension() {
        return 1;
    }

    @Override
    public double[] getGradientLogDensity(Object x) {
        double input = (double) x;
        if (input > lowerBound && input < upperBound) {
            return new double[]{0.0};
        } else if (input < lowerBound) {
            final double thetaL = (1.0 - lowerProbability - upperProbability)*lowerBound/(lowerProbability *(upperBound-lowerBound));
            return new double[]{(thetaL - 1) / input};
        } else {
            final double thetaR = (1.0 - lowerProbability - upperProbability) / (upperProbability * (upperBound - lowerBound));
            return new double[]{-thetaR};
        }
    }
}
