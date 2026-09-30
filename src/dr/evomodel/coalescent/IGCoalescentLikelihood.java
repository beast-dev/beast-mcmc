/*
 * IGCoalescentLikelihood.java
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

package dr.evomodel.coalescent;

import dr.evolution.coalescent.IntervalList;
import dr.evolution.coalescent.IntervalType;
import dr.evolution.util.Units;
import dr.evomodelxml.coalescent.CoalescentLikelihoodParser;
import dr.math.GammaFunction;

/**
 * A likelihood function for the coalescent that places an inverse-Gamma(alpha, beta) prior on
 * the (constant) population size and integrates it out analytically.
 *
 * The intervals of a genealogy are conditionally independent given the population size, but
 * become dependent once the population size is marginalised away, since they all share it.
 * The joint marginal density therefore cannot be obtained by multiplying per-interval marginal
 * densities; it must be computed from the sufficient statistics m (number of coalescent events)
 * and A (the sum of kChoose2 * duration over every interval, coalescent or not).
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @author Luiz Max Carvalho
 */
public final class IGCoalescentLikelihood extends AbstractCoalescentLikelihood implements Units {

    // PUBLIC STUFF
    public IGCoalescentLikelihood(IntervalList intervalList, double alpha, double beta) {

        super(CoalescentLikelihoodParser.COALESCENT_LIKELIHOOD, intervalList);

        this.alpha = alpha;
        this.beta = beta;

        this.coalescentEventStatisticValues = new double[getNumberOfCoalescentEvents()];
    }

    /**
     * Calculates the log likelihood of this set of coalescent intervals,
     * given doubles alpha and beta
     */
    protected double calculateLogLikelihood() {
        return calculateLogLikelihood(getIntervalList(), alpha, beta);
    }

    /**
     * Joint marginal log-likelihood of the intervals under a constant population size theta
     * given an IG(alpha, beta) prior on theta, with theta integrated out.
     */
    public static double calculateLogLikelihood(IntervalList intervals, double alpha, double beta) {

        double A = 0.0;
        int m = 0;
        double logCoalescentCoefficients = 0.0;

        final int n = intervals.getIntervalCount();
        for (int i = 0; i < n; i++) {

            final double duration = intervals.getInterval(i);
            final int lineageCount = intervals.getLineageCount(i);
            final double kChoose2 = lineageCount * (lineageCount - 1.0) / 2.0;

            A += kChoose2 * duration;

            if (intervals.getIntervalType(i) == IntervalType.COALESCENT) {
                m += 1;
                logCoalescentCoefficients += Math.log(kChoose2);
            }
        }

        return logCoalescentCoefficients
                + alpha * Math.log(beta)
                - GammaFunction.lnGamma(alpha)
                + GammaFunction.lnGamma(alpha + m)
                - (alpha + m) * Math.log(beta + A);
    }

    /**
     * Marginal log-density of a single coalescent interval of duration t among "coeff" = kChoose2
     * lineage pairs, with the population size integrated out under an IG(alpha, beta) prior. This
     * is correct in isolation (a Lomax density), but intervals share theta, so these must not be
     * summed across a tree in place of {@link #calculateLogLikelihood(IntervalList, double, double)}.
     */
    public static double marginalLogDensity(double t, double coeff, double alpha, double beta) {
        double lnum = Math.log(coeff) + Math.log(alpha) + alpha * Math.log(beta);
        double ldenom = (alpha + 1) * Math.log(coeff * t + beta);
        return lnum - ldenom;
    }

    /**
     * Log-probability of no coalescence within a sampling interval of duration t among "coeff" =
     * kChoose2 lineage pairs, with the population size integrated out under an IG(alpha, beta)
     * prior. Correct in isolation; see the caveat on {@link #marginalLogDensity}.
     */
    public static double invGammaNoCoalescentlogProb(double x, double coeff, double alpha, double beta) {
        return alpha * (Math.log(beta) - Math.log(coeff * x + beta));
    }

    // **************************************************************
    // CoalescentIntervalProvider IMPLEMENTATION
    // **************************************************************

    @Override
    public int getNumberOfCoalescentEvents() {
        return getIntervalList().getIntervalCount() / 2;
    }

    @Override
    public double getCoalescentEventsStatisticValue(int i) {
        if (i == 0) {
            IntervalList intervals = getIntervalList();
            final int intervalCount = intervals.getIntervalCount();
            for (int j = 0; j < coalescentEventStatisticValues.length; j++) {
                coalescentEventStatisticValues[j] = 0.0;
            }
            int counter = 0;
            for (int j = 0; j < intervalCount; j++) {
                final double kChoose2 = intervals.getLineageCount(j) * (intervals.getLineageCount(j) - 1.0) / 2.0;
                coalescentEventStatisticValues[counter] += intervals.getInterval(j) * kChoose2;
                if (intervals.getIntervalType(j) == IntervalType.COALESCENT) {
                    counter++;
                }
            }
        }
        return coalescentEventStatisticValues[i];
    }

    @Override
    public Type getUnits() {
        // TODO Auto-generated method stub
        return null;
    }

    @Override
    public void setUnits(Type units) {
        // TODO Auto-generated method stub
    }

    // PRIVATE STUFF
    private final double alpha;
    private final double beta;
    private final double[] coalescentEventStatisticValues;
}
