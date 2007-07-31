/*
 * FrequencyModel.java
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

package dr.evomodel.substmodel;

import dr.evolution.alignment.PatternList;
import dr.evolution.datatype.*;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.xml.*;
import dr.evoxml.DataTypeUtils;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

import java.util.logging.Logger;

/**
 * A model of equlibrium frequencies
 *
 * @version $Id: FrequencyModel.java,v 1.26 2005/05/24 20:25:58 rambaut Exp $
 *
 * @author Alexei Drummond
 * @author Andrew Rambaut
 */
public class FrequencyModel extends AbstractModel {

    public static final String FREQUENCIES = "frequencies";
    public static final String FREQUENCY_MODEL = "frequencyModel";

    public FrequencyModel(DataType dataType, Parameter frequencyParameter) {

        super(FREQUENCY_MODEL);
        this.frequencyParameter = frequencyParameter;
        addParameter(frequencyParameter);
        frequencyParameter.addBounds(new Parameter.DefaultBounds(1.0, 0.0, dataType.getStateCount()));
        this.dataType = dataType;
    }

    public final void setFrequency(int i, double value) {
        frequencyParameter.setParameterValue(i, value);
    }

    public final double getFrequency(int i) { return frequencyParameter.getParameterValue(i); }

    public int getFrequencyCount() { return frequencyParameter.getDimension(); }

    public final double[] getFrequencies() {
        double[] frequencies = new double[getFrequencyCount()];
        for (int i =0; i < frequencies.length; i++) {
            frequencies[i] = getFrequency(i);
        }
        return frequencies;
    }

    public DataType getDataType() { return dataType; }

    // *****************************************************************
    // Interface Model
    // *****************************************************************

    protected void handleModelChangedEvent(Model model, Object object, int index) {
        // no intermediates need recalculating....
    }

    protected void handleParameterChangedEvent(Parameter parameter, int index) {
        // no intermediates need recalculating....
    }

    protected void storeState() {} // no state apart from parameters to store
    protected void restoreState() {} // no state apart from parameters to restore
    protected void acceptState() {} // no state apart from parameters to accept
    protected void adoptState(Model source) {} // no state apart from parameters to adopt

    public Element createElement(Document doc) {
        throw new RuntimeException("Not implemented!");
    }

    /**
     * Reads a frequency model from an XMLObject.
     */
    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() { return FREQUENCY_MODEL; }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            DataType dataType = DataTypeUtils.getDataType(xo);

            Parameter freqsParam = (Parameter)xo.getSocketChild(FREQUENCIES);
            double[] frequencies = null;

            for (int i =0; i < xo.getChildCount(); i++) {
                Object obj = xo.getChild(i);
                if (obj instanceof PatternList) {
                    frequencies = ((PatternList)obj).getStateFrequencies();
                }
            }

            StringBuffer sb = new StringBuffer("Creating state frequencies model: ");
            if (frequencies != null) {
                if (freqsParam.getDimension() != frequencies.length) {
                    throw new XMLParseException("dimension of frequency parameter and number of sequence states don't match!");
                }
                for (int j = 0; j < frequencies.length; j++) {
                    freqsParam.setParameterValue(j, frequencies[j]);
                }
                sb.append("Using emprical frequencies from data ");
            } else {
                sb.append("Initial frequencies ");
            }
            sb.append("= {");
            sb.append(freqsParam.getParameterValue(0));
            for (int j = 1; j < freqsParam.getDimension(); j++) {
                sb.append(", ");
                sb.append(freqsParam.getParameterValue(j));
            }
            sb.append("}");
            Logger.getLogger("dr.evomodel").info(sb.toString());

            return new FrequencyModel(dataType, freqsParam);
        }

        public String getParserDescription() {
            return "A model of equilibrium base frequencies.";
        }

        public Class getReturnType() { return FrequencyModel.class; }

        public XMLSyntaxRule[] getSyntaxRules() { return rules; }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
                new XORRule(
                        new StringAttributeRule(DataType.DATA_TYPE, "The type of sequence data", DataType.getRegisteredDataTypeNames(), false),
                        new ElementRule(DataType.class)
                ),
                new ElementRule(FREQUENCIES,
                        new XMLSyntaxRule[] { new ElementRule(Parameter.class) })

        };
    };

    private DataType dataType = null;
    private Parameter frequencyParameter = null;

}
