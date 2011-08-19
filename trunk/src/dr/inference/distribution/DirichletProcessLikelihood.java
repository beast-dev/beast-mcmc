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
import dr.inference.model.Variable;
import dr.math.Binomial;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

/**
 * A class that returns the log likelihood of a set of discrete entities
 * being distributed in k classes according to a Dirichlet process.
 *
 * @author Andrew Rambaut
 * @author Trevor Bedford
 * @version $Id: BinomialLikelihood.java,v 1.5 2005/05/24 20:25:59 rambaut Exp $
 */

public class DirichletProcessLikelihood extends AbstractModelLikelihood {

    public static final String DIRICHLET_PROCESS_LIKELIHOOD = "dirichletProcessLikelihood";

    public DirichletProcessLikelihood(Parameter etaParameter, Parameter chiParameter) {

        super(DIRICHLET_PROCESS_LIKELIHOOD);

        this.etaParameter = etaParameter;
        this.chiParameter = chiParameter;
        addVariable(etaParameter);
        addVariable(chiParameter);

        K = etaParameter.getDimension();
        int count = 0;
        for (int i = 0; i < K; i++) {
            count += (int)etaParameter.getParameterValue(i);
        }
        N = count;
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
        double chi = chiParameter.getParameterValue(0);

        double logEtaj = 0;
        for (int j = 0; j < K; j++) {
            int eta = (int)etaParameter.getParameterValue(j) - 1;
            double logFactorial = 0;
            for (int k = 1; k <= eta; k++) {
                logFactorial += Math.log(k);
            }
            logEtaj += logFactorial;
        }
        double logDenominator = 0;
        for (int i = 1; i <= N; i++) {
            logDenominator += Math.log(chi + i - 1);
        }

        double logP = K * Math.log(chi) + logEtaj - logDenominator;

        return logP;
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

    private final Parameter etaParameter;
    private final Parameter chiParameter;
    private final int N, K;
}

