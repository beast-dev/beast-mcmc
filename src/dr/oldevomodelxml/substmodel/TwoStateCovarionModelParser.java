/*
 * TwoStateCovarionModelParser.java
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

import dr.evolution.datatype.TwoStateCovarion;
import dr.oldevomodel.substmodel.AbstractCovarionDNAModel;
import dr.oldevomodel.substmodel.FrequencyModel;
import dr.oldevomodel.substmodel.TwoStateCovarionModel;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * Parses an element from an DOM document into a TwoStateCovarionModel
 */
public class TwoStateCovarionModelParser extends AbstractXMLObjectParser {
    public static final String COVARION_MODEL = "covarionModel";
    public static final String ALPHA = "alpha";
    public static final String SWITCHING_RATE = "switchingRate";

    public String getParserName() {
        return COVARION_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Parameter alphaParameter;
        Parameter switchingRateParameter;

        XMLObject cxo = xo.getChild(AbstractCovarionDNAModel.FREQUENCIES);
        FrequencyModel freqModel = (FrequencyModel) cxo.getChild(FrequencyModel.class);

        TwoStateCovarion dataType = TwoStateCovarion.INSTANCE;  // fancy new datatype courtesy of Helen

        cxo = xo.getChild(ALPHA);
        alphaParameter = (Parameter) cxo.getChild(Parameter.class);

        // alpha must be positive and less than 1.0 because the fast rate is normalized to 1.0
        alphaParameter.addBounds(new Parameter.DefaultBounds(1.0, 0.0, 1));

        cxo = xo.getChild(SWITCHING_RATE);
        switchingRateParameter = (Parameter) cxo.getChild(Parameter.class);

        if (dataType != freqModel.getDataType()) {
            throw new XMLParseException("Data type of " + getParserName() + " element does not match that of its frequencyModel.");
        }

        TwoStateCovarionModel model = new TwoStateCovarionModel(dataType, freqModel, alphaParameter, switchingRateParameter);

        System.out.println(model);

        return model;
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "A covarion substitution model on binary data and a hidden rate state with two rates.";
    }

    public Class getReturnType() {
        return TwoStateCovarionModel.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new ElementRule(AbstractCovarionDNAModel.FREQUENCIES, FrequencyModel.class),
            new ElementRule(ALPHA, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class, true)}
            ),
            new ElementRule(SWITCHING_RATE, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class, true)}
            ),
    };
}
