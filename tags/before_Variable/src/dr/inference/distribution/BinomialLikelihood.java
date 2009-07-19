/*
 * BinomialLikelihood.java
 *
 * Copyright (C) 2002-2009 Alexei Drummond and Andrew Rambaut
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
 * BEAST is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
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

    public BinomialLikelihood(Parameter trialsParameter, Parameter proportionParameter, int[] counts) {

        super(BINOMIAL_LIKELIHOOD);

        this.trialsParameter = trialsParameter;
        this.proportionParameter = proportionParameter;
        addParameter(trialsParameter);
        addParameter(proportionParameter);
        this.counts = counts;

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

        double p = proportionParameter.getParameterValue(0);
        if (p <= 0 || p >= 1) return Double.NEGATIVE_INFINITY;

        double logP = Math.log(p);
        double log1MinusP = Math.log(1.0 - p);

        double logL = 0.0;
        for (int i = 0; i < trialsParameter.getDimension(); i++) {
            int trials = (int) Math.round(trialsParameter.getParameterValue(i));

            if (counts[i] > trials) return Double.NEGATIVE_INFINITY;
            logL += binomialLogLikelihood(trials, counts[i], logP, log1MinusP);
        }

        return logL;
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

    protected final void handleParameterChangedEvent(Parameter parameter, int index, Parameter.ChangeType type) {
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

    Binomial binom = new Binomial();
    Parameter trialsParameter;
    Parameter proportionParameter;
    int[] counts;
}

