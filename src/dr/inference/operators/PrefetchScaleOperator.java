/*
 * ScaleOperator.java
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

package dr.inference.operators;

import dr.inference.model.Bounds;
import dr.inference.model.Parameter;
import dr.inference.model.PrefetchableLikelihood;
import dr.inferencexml.operators.ScaleOperatorParser;
import dr.math.MathUtils;

/**
 * A generic scale operator for use with a multi-dimensional parameters.
 * Either scale all dimentions at once or scale one dimention at a time.
 * An optional bit vector and a threshold is used to vary the rate of the individual dimentions according
 * to their on/off status. For example a threshold of 1 means pick only "on" dimentions.
 *
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @version $Id: ScaleOperator.java,v 1.20 2005/06/14 10:40:34 rambaut Exp $
 */
public class PrefetchScaleOperator extends AbstractCoercableOperator implements PrefetchableOperator {

    private static final boolean PARALLEL_PREFETCH = false;
    private static final boolean PREFETCH_DEBUG = false;

    private Parameter parameter = null;
    private double scaleFactor = 0.5;

    private final PrefetchableLikelihood prefetchableLikelihood;
    private final int prefetchCount;
    private int currentPrefetch;
    private final double[] hastingsRatios;
    private final double[] draws;

    public PrefetchScaleOperator(PrefetchableLikelihood prefetchableLikelihood, Parameter parameter, double scaleFactor, CoercionMode mode, double weight) {

        super(mode);

        setWeight(weight);


        this.parameter = parameter;
        this.scaleFactor = scaleFactor;

        this.prefetchableLikelihood = prefetchableLikelihood;
        this.prefetchCount = prefetchableLikelihood.getPrefetchCount();
        hastingsRatios = new double[prefetchCount];
        draws = new double[prefetchCount];

        currentPrefetch = -1;

        if (parameter.getSize() > 1) {
            throw new IllegalArgumentException("PrefetchScaleOperator only works on parameters of dimension 1");
        }
    }

    @Override
    public double doOperation() {

        if (currentPrefetch < 0) {
            // Run the operator N times, cache the treeDataLikelihood evaluations.
            // store the hastings ratios. TreeDataLikelihood will also store likelihoods.
            currentPrefetch = 0;
            if (PREFETCH_DEBUG) {
                System.err.println("Prefetch: Drawing " + prefetchCount + " scaling factors");
            }

            Bounds<Double> bounds = parameter.getBounds();
            final double offset = bounds.getLowerLimit(0);

            // draw the prefetch operation scales
            for (int i = 0; i < prefetchCount; i++) {
                final double scale = (scaleFactor + (MathUtils.nextDouble() * ((1.0 / scaleFactor) - scaleFactor)));
                draws[i] = ((parameter.getParameterValue(0) - offset) * scale) + offset;
                hastingsRatios[i] = -Math.log(draws[i]);
            }

            if (PARALLEL_PREFETCH) {
                for (int i = 0; i < prefetchCount; i++) {
                    prefetchableLikelihood.startPrefetchOperation(i);

                    parameter.setParameterValue(0, draws[i]);
                }

                prefetchableLikelihood.prefetchLogLikelihoods();
            }
        }

        if (PREFETCH_DEBUG) {
            System.err.println("Prefetch: performing operation " + (currentPrefetch + 1));
        }

        if (PARALLEL_PREFETCH) {
            prefetchableLikelihood.startPrefetchOperation(currentPrefetch);
        } else {
            // A debugging option where there is no parallel processing of the operator
            // moves but they are simply done in sequence. The N operator instances
            // are drawn so that the random number sequence is conserved (allowing for
            // comparison with the parallel approach) and these are then just applied
            // in sequence as doOperation is called.

            parameter.setParameterValueQuietly(0, draws[currentPrefetch]);
        }

        return hastingsRatios[currentPrefetch];
    }


    @Override
    public void accept(double deviation) {
        super.accept(deviation);

        if (PARALLEL_PREFETCH) {
            prefetchableLikelihood.acceptPrefetch(currentPrefetch);
            parameter.setValue(0, draws[currentPrefetch]);
        }

        if (PREFETCH_DEBUG) {
            System.err.println("Prefetch: accepted");
        }

        currentPrefetch = -1;
    }

    @Override
    public void reject() {
        super.reject();

        // move on to the next prefetched
        currentPrefetch += 1;
        if (currentPrefetch == prefetchCount) {
            currentPrefetch = -1;
        }

        if (PREFETCH_DEBUG) {
            System.err.println("Prefetch: rejected");
        }

    }

    /**
     * This method should be overridden by operators that need to do something just before the return of doOperation.
     *
     * @param newValue the proposed parameter value
     * @param oldValue the old parameter value
     */
    void cleanupOperation(double newValue, double oldValue) {
        // DO NOTHING
    }

    public final String getOperatorName() {
        return "prefetchScale(" + parameter.getParameterName() + ")";
    }

    public double getCoercableParameter() {
        return Math.log(1.0 / scaleFactor - 1.0);
    }

    public void setCoercableParameter(double value) {
        scaleFactor = 1.0 / (Math.exp(value) + 1.0);
    }

    public double getRawParameter() {
        return scaleFactor;
    }

    public double getScaleFactor() {
        return scaleFactor;
    }

    public double getTargetAcceptanceProbability() {
        return 0.234;
    }

    public final String getPerformanceSuggestion() {

        double prob = Utils.getAcceptanceProbability(this);
        double targetProb = getTargetAcceptanceProbability();
        dr.util.NumberFormatter formatter = new dr.util.NumberFormatter(5);
        double sf = OperatorUtils.optimizeScaleFactor(scaleFactor, prob, targetProb);
        if (prob < getMinimumGoodAcceptanceLevel()) {
            return "Try setting scaleFactor to about " + formatter.format(sf);
        } else if (prob > getMaximumGoodAcceptanceLevel()) {
            return "Try setting scaleFactor to about " + formatter.format(sf);
        } else return "";
    }

    public String toString() {
        return ScaleOperatorParser.SCALE_OPERATOR + "(" + parameter.getParameterName() + " [" + scaleFactor + ", " + (1.0 / scaleFactor) + "]";
    }

    //PrefetchableOperator INTERFACE

    @Override
    public boolean prefetchingDone() {
        return currentPrefetch == prefetchCount - 1;
    }

    //PRIVATE STUFF

}
