/*
 * LogRateSubstitutionModelParser.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
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
 *
 */

package dr.evomodelxml.substmodel;

import dr.evolution.datatype.DataType;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.GlmSubstitutionModel;
import dr.evomodel.substmodel.LogAdditiveCtmcRateProvider;
import dr.evomodel.substmodel.LogRateSubstitutionModel;
import dr.evoxml.util.DataTypeUtils;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * @author Filippo Monti
 * @author Marc A. Suchard
 */

public class LogRateSubstitutionModelParser extends AbstractXMLObjectParser {

    public static final String LOG_RATE_SUBSTITUTION_MODEL = "logRateSubstitutionModel";
    private static final String NORMALIZE = "normalize";
    private static final String LOG_RATES = "logRates";

    public String getParserName() {
        return LOG_RATE_SUBSTITUTION_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        DataType dataType = DataTypeUtils.getDataType(xo);

        if (dataType == null) dataType = (DataType) xo.getChild(DataType.class);

        int rateCount = (dataType.getStateCount() - 1) * dataType.getStateCount(); // TODO not used anymore

        Parameter logRates = (Parameter) xo.getElementFirstChild(LOG_RATES);

        int length = logRates.getDimension();
        if (length != rateCount) {
            throw new XMLParseException("Rates parameter in " + getParserName() + " element should have " + (rateCount) + " dimensions. However the log rates dimension is " + length + ".");
        }

        LogAdditiveCtmcRateProvider lrm = new LogAdditiveCtmcRateProvider.DataAugmented.Basic(logRates.getId(), logRates);

        XMLObject cxo = xo.getChild(dr.oldevomodelxml.substmodel.ComplexSubstitutionModelParser.ROOT_FREQUENCIES);
        FrequencyModel rootFreq = (FrequencyModel) cxo.getChild(FrequencyModel.class);

        if (dataType != rootFreq.getDataType()) {
            throw new XMLParseException("Data type of " + getParserName() + " element does not match that of its rootFrequencyModel.");
        }

        boolean normalize = xo.getAttribute(NORMALIZE, true);

        LogRateSubstitutionModel model = new LogRateSubstitutionModel(xo.getId(), dataType, rootFreq, lrm);
        model.setNormalization(normalize);

        return model;
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "A general model of sequence substitution for any data type where the rates come from the generalized linear model.";
    }

    public Class getReturnType() {
        return GlmSubstitutionModel.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new XORRule(
                    new StringAttributeRule(DataType.DATA_TYPE, "The type of sequence data",
                            DataType.getRegisteredDataTypeNames(), false),
                    new ElementRule(DataType.class)
            ),
            new ElementRule(ComplexSubstitutionModelParser.ROOT_FREQUENCIES, FrequencyModel.class),
            new ElementRule(LOG_RATES, Parameter.class),
            AttributeRule.newBooleanRule(NORMALIZE, true),
    };
}
