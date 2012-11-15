/*
 * GeneralSubstitutionModel.java
 *
 * Copyright (C) 2002-2012 Alexei Drummond, Andrew Rambaut & Marc A. Suchard
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

package dr.app.beagle.evomodel.substmodel;

import dr.evolution.datatype.DataType;
import dr.inference.model.Parameter;
import dr.inference.model.DuplicatedParameter;

/**
 * <b>A general model of sequence substitution</b>. A general reversible class for any
 * data type.
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @author Marc Suchard
 */
public class GeneralSubstitutionModel extends BaseSubstitutionModel {

    /**
     * the rate which the others are set relative to
     */
    protected int ratesRelativeTo;

    public GeneralSubstitutionModel(String name, DataType dataType, FrequencyModel freqModel,
                                    Parameter parameter, int relativeTo) {
        this(name, dataType, freqModel, parameter, relativeTo, null);

    }

    /**
     * constructor
     *
     * @param dataType the data type
     */
    public GeneralSubstitutionModel(String name, DataType dataType, FrequencyModel freqModel,
                                    Parameter parameter, int relativeTo, EigenSystem eigenSystem) {

        super(name, dataType, freqModel, eigenSystem);

        ratesParameter = parameter;
        if (ratesParameter != null) {
            addVariable(ratesParameter);
            if (!(ratesParameter instanceof DuplicatedParameter))
                ratesParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0,
                        ratesParameter.getDimension()));
        }
        setRatesRelativeTo(relativeTo);
    }

    /**
     * constructor
     *
     * @param dataType the data type
     */
    protected GeneralSubstitutionModel(String name, DataType dataType, FrequencyModel freqModel, int relativeTo) {

        super(name, dataType, freqModel,
                null);

        setRatesRelativeTo(relativeTo);
    }

    protected void frequenciesChanged() {
        // Nothing to precalculate
    }

    protected void ratesChanged() {
        // Nothing to precalculate
    }

    protected void setupRelativeRates(double[] rates) {
        for (int i = 0; i < rates.length; i++) {
            if (i == ratesRelativeTo) {
                rates[i] = 1.0;
            } else if (i < ratesRelativeTo) {
                rates[i] = ratesParameter.getParameterValue(i);
            } else {
                rates[i] = ratesParameter.getParameterValue(i - 1);
            }
        }
    }

    /**
     * set which rate the others are relative to
     */
    public void setRatesRelativeTo(int ratesRelativeTo) {
        this.ratesRelativeTo = ratesRelativeTo;
    }

    // *****************************************************************
    // Interface Model
    // *****************************************************************


    protected void storeState() {
    } // nothing to do

    /**
     * Restore the additional stored state
     */
    protected void restoreState() {
        updateMatrix = true;
    }

    protected void acceptState() {
    } // nothing to do

    /**
     * Parses an element from an DOM document into a DemographicModel. Recognises
     * ConstantPopulation and ExponentialGrowth.
     */

    protected Parameter ratesParameter = null;
}