/*
 * SVSComplexSubstitutionModel.java
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

package dr.evomodel.substmodel;

import dr.inference.model.*;
import dr.evolution.datatype.DataType;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Marc Suchard
 */
public class SVSComplexSubstitutionModel extends ComplexSubstitutionModel implements Likelihood,
        BayesianStochasticSearchVariableSelection {

    public SVSComplexSubstitutionModel(String name, DataType dataType, FrequencyModel freqModel,
                                       Parameter ratesParameter, Parameter indicatorsParameter) {
        super(name, dataType, freqModel, ratesParameter);

        if (indicatorsParameter == null) {
            this.indicatorsParameter = new Parameter.Default(ratesParameter.getDimension(), 1.0);
        } else {
            this.indicatorsParameter  = indicatorsParameter;
            addVariable(indicatorsParameter);
        }

        setupIndicatorDimensionNames(-1);
    }

    @Override
    protected void setupRelativeRates(double[] rates) {
        for (int i = 0; i < rates.length; i++) {
            rates[i] = ratesParameter.getParameterValue(i) * indicatorsParameter.getParameterValue(i);
        }
    }

    protected void setupIndicatorDimensionNames(int relativeTo) {
        List<String> indicatorNames = new ArrayList<String>();

        String indicatorPrefix = indicatorsParameter.getParameterName();

        for (int i = 0; i < dataType.getStateCount(); ++i) {
            for (int j = i + 1; j < dataType.getStateCount(); ++j) {
                indicatorNames.add(getDimensionString(i, j, indicatorPrefix));
            }
        }

        for (int j = 0; j < dataType.getStateCount(); ++j) {
            for (int i = j + 1; i < dataType.getStateCount(); ++i) {
                indicatorNames.add(getDimensionString(i, j, indicatorPrefix));
            }
        }

        String[] tmp = new String[0];
        indicatorsParameter.setDimensionNames(indicatorNames.toArray(tmp));
    }

    public Parameter getIndicators() {
        return indicatorsParameter;
    }

    public boolean validState() {
        return !updateMatrix ||
                BayesianStochasticSearchVariableSelection.Utils.connectedAndWellConditioned(probability,this);
    }

    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        if (variable == ratesParameter && indicatorsParameter.getParameterValue(index) == 0)
            return; // Does not affect likelihood
        super.handleVariableChangedEvent(variable,index,type);
    }

    public Model getModel() {
        return this;
    }

    /**
     * Forces a complete recalculation of the likelihood next time getLikelihood is called
     */
    public void makeDirty() {
        updateMatrix = true;
    }

    /**
     * @return A detailed name of likelihood for debugging.
     */
//    public String prettyName() {
//        return "SVSComplexSubstitutionModel-connectedness";
//    }

    @Override
    public boolean isUsed() {
        return super.isUsed() && isUsed;
    }

    public void setUsed() {
        isUsed = true;
    }

    private boolean isUsed = false;

    private double[] probability = null;

    private final Parameter indicatorsParameter;
}
