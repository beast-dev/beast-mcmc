/*
 * FrequencyModel.java
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
import dr.oldevomodelxml.substmodel.FrequencyModelParser;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

/**
 * A model of equlibrium frequencies
 *
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @author Marc Suchard
 * @version $Id: FrequencyModel.java,v 1.26 2005/05/24 20:25:58 rambaut Exp $
 */
public class FrequencyModel extends AbstractModel {
    /**
     * A constructor which allows a more programmatic approach with
     * fixed frequencies.
     * @param dataType              DataType
     * @param frequencyParameter    double[]
     */
    public FrequencyModel(DataType dataType, double[] frequencyParameter) {
        this(dataType, new Parameter.Default(frequencyParameter));
    }

    public FrequencyModel(DataType dataType, Parameter frequencyParameter) {

        super(FrequencyModelParser.FREQUENCY_MODEL);

        double sum = getSumOfFrequencies(frequencyParameter);

        if (Math.abs(sum - 1.0) > 1e-8) {
            throw new IllegalArgumentException("Frequencies do not sum to 1, they sum to " + sum);
        }

        this.frequencyParameter = frequencyParameter;
        addVariable(frequencyParameter);
        frequencyParameter.addBounds(new Parameter.DefaultBounds(1.0, 0.0, frequencyParameter.getDimension()));
        this.dataType = dataType;
    }

    /**
     * @param frequencies the frequencies
     * @return return the sum of frequencies
     */
    private double getSumOfFrequencies(Parameter frequencies) {
        double total = 0.0;
        for (int i = 0; i < frequencies.getDimension(); i++) {
            total += frequencies.getParameterValue(i);
        }
        return total;
    }

    public void setFrequency(int i, double value) {
        frequencyParameter.setParameterValue(i, value);
    }

    public double getFrequency(int i) {
        return frequencyParameter.getParameterValue(i);
    }

    public int getFrequencyCount() {
        return frequencyParameter.getDimension();
    }

    public Parameter getFrequencyParameter() {
        return frequencyParameter;
    }

    public double[] getFrequencies() {
        double[] frequencies = new double[getFrequencyCount()];
        for (int i = 0; i < frequencies.length; i++) {
            frequencies[i] = getFrequency(i);
        }
        return frequencies;
    }

    public double[] getCumulativeFrequencies() {
        double[] frequencies = getFrequencies();
        for (int i = 1; i < frequencies.length; i++) {
            frequencies[i] += frequencies[i - 1];
        }
        return frequencies;
    }

    public DataType getDataType() {
        return dataType;
    }

    // *****************************************************************
    // Interface Model
    // *****************************************************************

    protected void handleModelChangedEvent(Model model, Object object, int index) {
        // no intermediates need recalculating....
    }

    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        // no intermediates need recalculating....
    }

    protected void storeState() {
    } // no state apart from parameters to store

    protected void restoreState() {
    } // no state apart from parameters to restore

    protected void acceptState() {
    } // no state apart from parameters to accept

    public Element createElement(Document doc) {
        throw new RuntimeException("Not implemented!");
    }

    private DataType dataType = null;
    Parameter frequencyParameter = null;

}
