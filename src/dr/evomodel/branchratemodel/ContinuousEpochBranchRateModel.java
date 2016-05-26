/*
 * ContinuousEpochBranchRateModel.java
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

package dr.evomodel.branchratemodel;

import dr.inference.model.Parameter;
import dr.inference.model.Variable;

/**
 * @author Marc Suchard
 */
public class ContinuousEpochBranchRateModel extends RateEpochBranchRateModel {

    /**
     * The constructor. For an N-epoch model, there should be N rate paramters and N-1 transition times.
     *
     * @param timeParameters an array of transition time parameters
     * @param rateParameters an array of rate parameters
     */
    public ContinuousEpochBranchRateModel(Parameter[] timeParameters, Parameter[] rateParameters, Parameter rootHeight) {
        super(timeParameters, rateParameters);
        this.rootHeight = rootHeight;
        addVariable(rootHeight);
        normalizationKnown = false;
    }

    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        super.handleVariableChangedEvent(variable, index, type);
        normalizationKnown = false;
    }

    private void normalize() {
        normalization = 0.0;
        double startTime = 0.0;
        double endTime = rootHeight.getParameterValue(0);
        int j = 0;
        while( j < timeParameters.length && endTime > timeParameters[j].getParameterValue(0)) {
            final double nextTime = timeParameters[j].getParameterValue(0);
            normalization += (nextTime - startTime) * rateParameters[j].getParameterValue(0);
            startTime = nextTime;
            j++;
        }
        normalization += (endTime - startTime) * rateParameters[j].getParameterValue(0);
    }

    protected void storeState() {
        savedNormalization = normalization;
    }

    protected void restoreState() {
        normalization = savedNormalization;
    }

    protected double normalizeRate(double rate) {
        if (!normalizationKnown)
            normalize();
        return rate / normalization;
    }

    private Parameter rootHeight;
    private double normalization;
    private double savedNormalization;
    private boolean normalizationKnown = false;

}
