/*
 * BinaryCovarionModelParser.java
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

package dr.app.beagle.evomodel.parsers;

import dr.app.beagle.evomodel.substmodel.BinaryCovarionModel;
import dr.evolution.datatype.TwoStateCovarion;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * Parses an element from an DOM document into a TwoStateCovarionModel
 */
public class BinaryCovarionModelParser extends AbstractXMLObjectParser {
    public static final String COVARION_MODEL = dr.evomodelxml.substmodel.BinaryCovarionModelParser.COVARION_MODEL;
    public static final String ALPHA = dr.evomodelxml.substmodel.BinaryCovarionModelParser.ALPHA;
    public static final String SWITCHING_RATE = dr.evomodelxml.substmodel.BinaryCovarionModelParser.SWITCHING_RATE;
    public static final String FREQUENCIES = dr.evomodelxml.substmodel.BinaryCovarionModelParser.FREQUENCIES;
    public static final String HIDDEN_FREQUENCIES = dr.evomodelxml.substmodel.BinaryCovarionModelParser.HIDDEN_FREQUENCIES;
    public static final String VERSION = dr.evomodelxml.substmodel.BinaryCovarionModelParser.VERSION;
    public static final dr.evomodel.substmodel.BinaryCovarionModel.Version DEFAULT_VERSION =
            dr.evomodelxml.substmodel.BinaryCovarionModelParser.DEFAULT_VERSION;

    public String getParserName() {
        return COVARION_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Parameter alphaParameter;
        Parameter switchingRateParameter;

        XMLObject cxo = xo.getChild(FREQUENCIES);
        Parameter frequencies = (Parameter) cxo.getChild(Parameter.class);

        cxo = xo.getChild(HIDDEN_FREQUENCIES);
        Parameter hiddenFrequencies = (Parameter) cxo.getChild(Parameter.class);

        cxo = xo.getChild(ALPHA);
        alphaParameter = (Parameter) cxo.getChild(Parameter.class);

        // alpha must be positive and less than 1.0 because the fast rate is normalized to 1.0
        alphaParameter.addBounds(new Parameter.DefaultBounds(1.0, 0.0, 1));
        hiddenFrequencies.addBounds(new Parameter.DefaultBounds(1.0, 0.0, hiddenFrequencies.getDimension()));
        frequencies.addBounds(new Parameter.DefaultBounds(1.0, 0.0, frequencies.getDimension()));

        cxo = xo.getChild(SWITCHING_RATE);
        switchingRateParameter = (Parameter) cxo.getChild(Parameter.class);

        dr.evomodel.substmodel.BinaryCovarionModel.Version version = DEFAULT_VERSION;
        if (xo.hasAttribute(VERSION)) {
            version = dr.evomodel.substmodel.BinaryCovarionModel.Version.parseFromString(xo.getStringAttribute(VERSION));
        }

        BinaryCovarionModel model = new BinaryCovarionModel(TwoStateCovarion.INSTANCE,
                frequencies, hiddenFrequencies, alphaParameter, switchingRateParameter, version);

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
        return BinaryCovarionModel.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(FREQUENCIES, Parameter.class),
            new ElementRule(HIDDEN_FREQUENCIES, Parameter.class),
            new ElementRule(ALPHA,
                    new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class, true)}
            ),
            new ElementRule(SWITCHING_RATE,
                    new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class, true)}
            ),
            AttributeRule.newStringRule(VERSION, true),
    };


}
