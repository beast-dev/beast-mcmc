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
 * The truncated Cauchy distribution with soft bound in Inoue, Donoghue and Yang, 2010
 * reference: paml4.8 package by Dr. Ziheng Yang
 *
 * @author Xiang Ji
 * @author Marc Suchard
 */
public class TruncatedSoftBoundCauchyDistribution extends AbstractContinuousDistribution implements Distribution, GradientProvider, Citable {

    private final double lowerBound;
    private final double offSetProportion;
    private final double scaleMultiplier;
    private final double tailL;

    private final double A;
    private final double t0;
    private final double s;


    public TruncatedSoftBoundCauchyDistribution(double lowerBound,
                                                double offSetProportion,
                                                double scaleMultiplier) {
        this(lowerBound, offSetProportion, scaleMultiplier, 0.025);
    }

    public TruncatedSoftBoundCauchyDistribution(double lowerBound,
                                                double offSetProportion,
                                                double scaleMultiplier,
                                                double tailL){

        this.lowerBound = lowerBound;
        this.offSetProportion = offSetProportion;
        this.scaleMultiplier = scaleMultiplier;
        this.tailL = tailL;

        this.A = 0.5 + 1.0 / Math.PI * Math.atan2(scaleMultiplier, offSetProportion);
        this.t0 = lowerBound * (1.0 + offSetProportion);
        this.s = lowerBound * scaleMultiplier;
    }

    @Override
    public double pdf(double x) {
        assert x > 0;
        if (x > lowerBound) {
            final double z = (x - t0) / s;
            return (1-tailL)/(Math.PI*A*s*(1 + z*z));
        } else {
            final double z = offSetProportion / scaleMultiplier;
            final double thetaL =  (1/tailL-1) / (Math.PI*A*scaleMultiplier*(1 + z*z));
            return tailL*thetaL/lowerBound * Math.pow(x / lowerBound, thetaL - 1.0);
        }
    }

    /*
    Reference: line 2263 - line 2279 of mcmctree.c from paml4.8 package by Dr. Ziheng Yang
     */
    @Override
    public double logPdf(double x) {
        assert x > 0;
        if (x > lowerBound) {
            final double z = (x - t0) / s;
            return Math.log((1-tailL)/(Math.PI*A*s*(1 + z*z)));
        } else {
            final double z = offSetProportion / scaleMultiplier;
            final double thetaL =  (1.0 / tailL-1) / (Math.PI*A*scaleMultiplier*(1 + z*z));
            return Math.log(tailL * thetaL / lowerBound) + (thetaL-1)*Math.log(x / lowerBound);
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
        return "Truncated Cauchy distribution with soft bound.";
    }

    @Override
    public List<Citation> getCitations() {
        return Collections.singletonList(CITATION);
    }

    public static Citation CITATION = new Citation(
            new Author[]{
                    new Author("J", "Inoue"),
                    new Author("PCJ", "Donoghue"),
                    new Author("Z", "Yang")
            },
            "The impact of the representation of fossil calibrations on Bayesian estimation of species divergence times",
            2010,
            "Systematic Biology",
            59,
            74,
            89,
            Citation.Status.PUBLISHED
    );

    @Override
    public int getDimension() {
        return 1;
    }

    @Override
    public double[] getGradientLogDensity(Object x) {
        double input = (double) x;
        assert input > 0;
        if (input > lowerBound) {
            final double z = (input - t0) / s;
            return new double[]{-2.0 * z / (1 + z * z) / s};
        } else {
            final double z = offSetProportion / scaleMultiplier;
            final double thetaL =  (1.0 / tailL-1) / (Math.PI*A*scaleMultiplier*(1 + z*z));
            return new double[]{(thetaL - 1) / input};
        }
    }
}
