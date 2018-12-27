/*
 * MCMCOptions.java
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

package dr.inference.mcmc;

import dr.inference.markovchain.MarkovChain;

/**
 * A class that brings together the auxillary information associated
 * with an MCMC analysis.
 *
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @version $Id: MCMCOptions.java,v 1.7 2005/05/24 20:25:59 rambaut Exp $
 */
public class MCMCOptions {

    private final long chainLength;
    private final long fullEvaluationCount;
    private final int minOperatorCountForFullEvaluation;
    private final double evaluationTestThreshold;
    private final boolean useAdaptation;
    private final long adaptationDelay;
    private final double adaptationTarget;
    private final double temperature;

    /**
     * constructor
     * @param chainLength
     */
    public MCMCOptions(long chainLength) {
        this(chainLength, 2000, 1, MarkovChain.EVALUATION_TEST_THRESHOLD, true, 0, 0.234, 1.0);
    }

    /**
     * constructor
     * @param chainLength
     * @param fullEvaluationCount
     * @param minOperatorCountForFullEvaluation
     * @param evaluationTestThreshold
     * @param useAdaptation
     * @param adaptationDelay
     * @param temperature
     */
    public MCMCOptions(long chainLength, long fullEvaluationCount, int minOperatorCountForFullEvaluation,
                       double evaluationTestThreshold, boolean useAdaptation, long adaptationDelay, double adaptationTarget,
                       double temperature) {
        this.chainLength = chainLength;
        this.fullEvaluationCount = fullEvaluationCount;
        this.minOperatorCountForFullEvaluation = minOperatorCountForFullEvaluation;
        this.evaluationTestThreshold = evaluationTestThreshold;
        this.useAdaptation = useAdaptation;
        this.adaptationDelay = adaptationDelay;
        this.adaptationTarget = adaptationTarget;
        this.temperature = temperature;
    }

    /**
     * @return the chain length of the MCMC analysis
     */
    public final long getChainLength() {
        return chainLength;
    }

    public final long getFullEvaluationCount() {
        return fullEvaluationCount;
    }

    public double getEvaluationTestThreshold() {
        return evaluationTestThreshold;
    }

    public final boolean useAdaptation() {
        return useAdaptation;
    }

    public final long getAdaptationDelay() {
        return adaptationDelay;
    }

    public final double getAdaptationTarget() {
        return adaptationTarget;
    }

    public final double getTemperature() {
        return temperature;
    }

    public int minOperatorCountForFullEvaluation() {
        return minOperatorCountForFullEvaluation;
    }
}
