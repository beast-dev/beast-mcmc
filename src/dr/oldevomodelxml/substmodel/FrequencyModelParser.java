/*
 * FrequencyModelParser.java
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

package dr.oldevomodelxml.substmodel;

import dr.evolution.alignment.PatternList;
import dr.evolution.datatype.DataType;
import dr.oldevomodel.substmodel.FrequencyModel;
import dr.evoxml.util.DataTypeUtils;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.text.NumberFormat;
import java.util.logging.Logger;

/**
 * Reads a frequency model from an XMLObject.
 */
public class FrequencyModelParser extends AbstractXMLObjectParser {
    public static final String FREQUENCIES = "frequencies";
    public static final String FREQUENCY_MODEL = "frequencyModel";
    public static final String NORMALIZE = "normalize";


    public String[] getParserNames() {
        return new String[]{
                getParserName(), "beast_" + getParserName()
        };
    }

    public String getParserName() {
        return FREQUENCY_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        DataType dataType = DataTypeUtils.getDataType(xo);

        Parameter freqsParam = (Parameter) xo.getElementFirstChild(FREQUENCIES);
        double[] frequencies = null;

        for (int i = 0; i < xo.getChildCount(); i++) {
            Object obj = xo.getChild(i);
            if (obj instanceof PatternList) {
                frequencies = ((PatternList) obj).getStateFrequencies();
                break;
            }
        }

        StringBuilder sb = new StringBuilder("Creating state frequencies model '" + freqsParam.getParameterName() + "': ");
        if (frequencies != null) {
            if (freqsParam.getDimension() != frequencies.length) {
                throw new XMLParseException("dimension of frequency parameter and number of sequence states don't match!");
            }
            for (int j = 0; j < frequencies.length; j++) {
                freqsParam.setParameterValue(j, frequencies[j]);
            }
            sb.append("Using empirical frequencies from data ");
        } else {
            sb.append("Initial frequencies ");
        }
        sb.append("= {");

        double sum = 0;
        for (int j = 0; j < freqsParam.getDimension(); j++) {
            sum += freqsParam.getParameterValue(j);
        }

        if (xo.getAttribute(NORMALIZE, false)) {
            for (int j = 0; j < freqsParam.getDimension(); j++) {
                if (sum != 0)
                    freqsParam.setParameterValue(j, freqsParam.getParameterValue(j) / sum);
                else
                    freqsParam.setParameterValue(j, 1.0 / freqsParam.getDimension());
            }
            sum = 1.0;
        }

        if (Math.abs(sum - 1.0) > 1e-8) {
            throw new XMLParseException("Frequencies do not sum to 1 (they sum to " + sum + ")");
        }


        NumberFormat format = NumberFormat.getNumberInstance();
        format.setMaximumFractionDigits(5);

        sb.append(format.format(freqsParam.getParameterValue(0)));
        for (int j = 1; j < freqsParam.getDimension(); j++) {
            sb.append(", ");
            sb.append(format.format(freqsParam.getParameterValue(j)));
        }
        sb.append("}");
        Logger.getLogger("dr.evomodel").info(sb.toString());

        return new FrequencyModel(dataType, freqsParam);
    }

    public String getParserDescription() {
        return "A model of equilibrium base frequencies.";
    }

    public Class getReturnType() {
        return FrequencyModel.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newBooleanRule(NORMALIZE, true),

            new ElementRule(PatternList.class, "Initial value", 0, 1),

            new XORRule(
                    new StringAttributeRule(DataType.DATA_TYPE, "The type of sequence data",
                            DataType.getRegisteredDataTypeNames(), false),
                    new ElementRule(DataType.class)
            ),

            new ElementRule(FREQUENCIES,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),

    };


}
