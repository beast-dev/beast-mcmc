/*
 * GeneralSubstitutionModel.java
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

package dr.oldevomodel.substmodel;

import dr.evolution.datatype.DataType;
import dr.oldevomodelxml.substmodel.GeneralSubstitutionModelParser;
import dr.inference.model.DuplicatedParameter;
import dr.inference.model.Parameter;

/**
 * <b>A general model of sequence substitution</b>. A general reversible class for any
 * data type.
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: GeneralSubstitutionModel.java,v 1.37 2006/05/05 03:05:10 alexei Exp $
 */
public class GeneralSubstitutionModel extends AbstractSubstitutionModel implements dr.util.XHTMLable {

    /**
     * the rate which the others are set relative to
     */
    protected int ratesRelativeTo;

    /**
     * constructor
     *
     * @param dataType   the data type
     * @param freqModel  the equilibrium frequency model - this must match the data type
     * @param parameter  the rates parameter, minus the rate that they are specified relative to
     * @param relativeTo the index of the rate that all other are specified relative to
     */
    public GeneralSubstitutionModel(
            DataType dataType,
            FrequencyModel freqModel,
            Parameter parameter,
            int relativeTo) {

        super(GeneralSubstitutionModelParser.GENERAL_SUBSTITUTION_MODEL, dataType, freqModel);

        ratesParameter = parameter;
        if (ratesParameter != null) {
            addVariable(ratesParameter);
            if (!(ratesParameter instanceof DuplicatedParameter))
                ratesParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, ratesParameter.getDimension()));
        }
        setRatesRelativeTo(relativeTo);
    }

    /**
     * constructor
     *
     * @param name       a name to give the substitution model
     * @param dataType   the data type
     * @param freqModel  the equilibrium frequency model - this must match the data type
     * @param relativeTo the index of the rate that all other are specified relative to
     */
    protected GeneralSubstitutionModel(
            String name,
            DataType dataType,
            FrequencyModel freqModel,
            int relativeTo) {

        super(name, dataType, freqModel);

        setRatesRelativeTo(relativeTo);
    }

    protected void frequenciesChanged() {
        // Nothing to precalculate
    }

    protected void ratesChanged() {
        // Nothing to precalculate
    }

    protected void setupRelativeRates() {

        for (int i = 0; i < relativeRates.length; i++) {
            if (i == ratesRelativeTo) {
                relativeRates[i] = 1.0;
            } else if (i < ratesRelativeTo) {
                relativeRates[i] = ratesParameter.getParameterValue(i);
            } else {
                relativeRates[i] = ratesParameter.getParameterValue(i - 1);
            }
        }
    }

    /**
     * set which rate the others are relative to
     *
     * @param ratesRelativeTo the index of the rate in the matrix that all other
     *                        rates are parameterized relative to.
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

    // **************************************************************
    // XHTMLable IMPLEMENTATION
    // **************************************************************

    public String toXHTML() {
        StringBuffer buffer = new StringBuffer();

        buffer.append("<em>General Model</em>");

        return buffer.toString();
    }

    protected Parameter ratesParameter = null;
}
