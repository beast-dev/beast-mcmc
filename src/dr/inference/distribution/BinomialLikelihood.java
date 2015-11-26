/*
 * BinomialLikelihood.java
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

import dr.inference.model.AbstractModelLikelihood;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.math.Binomial;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

/**
 * A class that returns the log likelihood of a set of data (statistics)
 * being distributed according to a binomial distribution.
 *
 * @author Alexei Drummond
 * @version $Id: BinomialLikelihood.java,v 1.5 2005/05/24 20:25:59 rambaut Exp $
 */

public class BinomialLikelihood extends AbstractModelLikelihood {

    public static final String BINOMIAL_LIKELIHOOD = "binomialLikelihood";

    public BinomialLikelihood(Parameter trialsParameter, Parameter proportionParameter, Parameter countsParameter,
                              boolean onLogitScale) {

        super(BINOMIAL_LIKELIHOOD);

        this.trialsParameter = trialsParameter;
        this.proportionParameter = proportionParameter;
        this.countsParameter = countsParameter;
        addVariable(trialsParameter);
        addVariable(proportionParameter);
        addVariable(countsParameter);

        this.onLogitScale = onLogitScale;

    }

    // **************************************************************
    // Likelihood IMPLEMENTATION
    // **************************************************************


    public Model getModel() {
        return this;
    }

    /**
     * Calculate the log likelihood of the current state.
     *
     * @return the log likelihood.
     */
    public double getLogLikelihood() {

        // Get first success probability statistics
        SufficientStatistics ss = new SufficientStatistics(proportionParameter.getParameterValue(0), onLogitScale);
        if (ss.outOfBounds) return Double.NEGATIVE_INFINITY;

        final boolean hasMultipleProportions = proportionParameter.getDimension() > 1;

        double logL = 0.0;
        for (int i = 0; i < trialsParameter.getDimension(); i++) {

            if (hasMultipleProportions && i > 0) {
                ss = new SufficientStatistics(proportionParameter.getParameterValue(i), onLogitScale);
                if (ss.outOfBounds) return Double.NEGATIVE_INFINITY;
            }

            int trials = (int) Math.round(trialsParameter.getParameterValue(i));
            int counts = (int) Math.round(countsParameter.getParameterValue(i));

            if (counts > trials) return Double.NEGATIVE_INFINITY;
            logL += binomialLogLikelihood(trials, counts, ss.logP, ss.log1MinusP);
        }
        return logL;
    }

    class SufficientStatistics {

        double logP;
        double log1MinusP;
        boolean outOfBounds;

        SufficientStatistics(double theta, boolean onLogitScale) {
            if (onLogitScale) {
                // e(x) / (1 + e(x)) and 1  / (1 + e(x))
                final double logDenominator = Math.log(1.0 + Math.exp(theta));
                logP = theta - logDenominator;
                log1MinusP = -logDenominator;
                outOfBounds = false;
            } else {
                logP = Math.log(theta);
                log1MinusP = Math.log(1.0 - theta);
                outOfBounds = (theta <= 0 || theta >= 1);
            }
        }
    }

    public void makeDirty() {
    }

    public void acceptState() {
        // DO NOTHING
    }

    public void restoreState() {
        // DO NOTHING
    }

    public void storeState() {
        // DO NOTHING
    }

    protected void handleModelChangedEvent(Model model, Object object, int index) {
        // DO NOTHING
    }

    protected final void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        // DO NOTHING
    }

    /**
     * @return the binomial likelihood of obtaining the gicen count in the given number of trials,
     *         when the log of the probability is logP.
     */
    private double binomialLogLikelihood(int trials, int count, double logP, double log1MinusP) {
        return Math.log(Binomial.choose(trials, count)) + (logP * count) + (log1MinusP * (trials - count));
    }

    // **************************************************************
    // XMLElement IMPLEMENTATION
    // **************************************************************

    public Element createElement(Document d) {
        throw new RuntimeException("Not implemented yet!");
    }

    Parameter trialsParameter;
    Parameter proportionParameter;
    Parameter countsParameter;

    private final boolean onLogitScale;
}

