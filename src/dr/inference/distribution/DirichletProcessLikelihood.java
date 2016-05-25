/*
 * DirichletProcessLikelihood.java
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

import dr.inference.model.*;

/**
 * A class that returns the log likelihood of a set of discrete entities
 * being distributed in k classes according to a Dirichlet process.
 *
 * @author Andrew Rambaut
 * @author Trevor Bedford
 * @author Marc Suchard
 * @version $Id: DirichletProcessLikelihood.java,v 1.7 2005/05/24 20:25:59 rambaut Exp $
 */

public class DirichletProcessLikelihood extends AbstractModelLikelihood {

    public static final String DIRICHLET_PROCESS_LIKELIHOOD = "dirichletProcessLikelihood";

    public DirichletProcessLikelihood(Statistic etaParameter, Parameter chiParameter) {

        super(DIRICHLET_PROCESS_LIKELIHOOD);

        this.etaParameter = etaParameter;
        this.chiParameter = chiParameter;
        addVariable(chiParameter);

        K = etaParameter.getDimension();
        int count = 0;
        for (int i = 0; i < K; i++) {
            count += (int)etaParameter.getStatisticValue(i);
        }
        N = count;

        // create a look up table for all log factorials up to N
        logFactorials = new double[N];
        logFactorials[0] = 0.0;
        for (int j = 1; j < N; j++) {
//            logFactorials[j] = 0; // the log factorial for 0
//            for (int k = 1; k <= j; k++) {
//                logFactorials[j] += Math.log(k);
//            }
            logFactorials[j] = logFactorials[j - 1] + Math.log(j);
        }
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
        int K1 = 0;
        for (int j = 0; j < K; j++) {
            int eta = (int)etaParameter.getStatisticValue(j);
            if (eta > N) {
                throw new RuntimeException("Illegal eta value");
            }
            if (eta > 0) {
//                double logFactorial = 0;
//                for (int k = 1; k <= (eta - 1); k++) {
//                    logFactorial += Math.log(k);
//                }
                double logFactorial = logFactorials[eta - 1];
                logEtaj += logFactorial;
                // count the number of actually occupied classes
                K1++;
            }
        }
        double logDenominator = 0;
        for (int i = 1; i <= N; i++) {
            logDenominator += Math.log(chi + i - 1);
        }

        double logP = K1 * Math.log(chi) + logEtaj - logDenominator;

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

    public Statistic getEtaParameter() {
        return etaParameter;
    }

    public Parameter getChiParameter() {
        return chiParameter;
    }

    public int getN() {
        return N;
    }

    private final Statistic etaParameter;
    private final Parameter chiParameter;
    private final int N, K;
    private final double[] logFactorials;
}

