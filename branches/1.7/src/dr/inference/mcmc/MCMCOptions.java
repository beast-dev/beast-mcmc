/*
 * MCMCOptions.java
 *
 * Copyright (C) 2002-2006 Alexei Drummond and Andrew Rambaut
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

/**
 * A class that brings together the auxillary information associated
 * with an MCMC analysis.
 *
 * @author Alexei Drummond
 * @version $Id: MCMCOptions.java,v 1.7 2005/05/24 20:25:59 rambaut Exp $
 */
public class MCMCOptions {

    private long chainLength;
    private int fullEvaluationCount = 2000;
    private int minOperatorCountForFullEvaluation = 1;
    private boolean coercion = true;
    private int coercionDelay = 0;
    private double temperature = 1.0;

    public MCMCOptions() {
    }

    /**
     * @return the chain length of the MCMC analysis
     */
    public final long getChainLength() {
        return chainLength;
    }

    public final int fullEvaluationCount() {
        return fullEvaluationCount;
    }

    public final boolean useCoercion() {
        return coercion;
    }


    public final int getCoercionDelay() {
        return coercionDelay;
    }

    public final double getTemperature() {
        return temperature;
    }

    public final void setTemperature(double temperature) {
        this.temperature = temperature;
    }

    public final void setChainLength(long length) {
        chainLength = length;
    }

    public final void setFullEvaluationCount(int fullEvaluationCount) {
        this.fullEvaluationCount = fullEvaluationCount;
    }

    public final void setUseCoercion(boolean coercion) {
        this.coercion = coercion;
        if (!coercion) coercionDelay = 0;
    }

    public final void setCoercionDelay(int coercionDelay) {
        this.coercionDelay = coercionDelay;
    }

    public int minOperatorCountForFullEvaluation() {
        return minOperatorCountForFullEvaluation;
    }

    public final void setMinOperatorCountForFullEvaluation(int count) {
        this.minOperatorCountForFullEvaluation = count;
    }
}
