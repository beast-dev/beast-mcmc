/*
 * MicrosatelliteModelSelectOperatorParser.java
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

package dr.evomodelxml.operators;

import dr.inference.model.Parameter;
import dr.inference.operators.MCMCOperator;
import dr.evomodel.operators.MicrosatelliteModelSelectOperator;
import dr.xml.*;

/**
 * Parser for MicrosatelliteModelSelectOperatorParser
 */
public class MicrosatelliteModelSelectOperatorParser extends AbstractXMLObjectParser {

    public static final String MODEL_INDICATORS = "modelIndicators";
    public static final String MODEL_CHOOSE = "modelChoose";

    public String getParserName() {
        return "msatModelSelectOperator";
    }
         public Object parseXMLObject(XMLObject xo) throws XMLParseException {
             double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);
        Parameter modelChoose = (Parameter)xo.getElementFirstChild(MODEL_CHOOSE);
        XMLObject xoInd = xo.getChild(MODEL_INDICATORS);
             int childNum = xoInd.getChildCount();
        System.out.println("There are 12 potential models");
        Parameter[] modelIndicators = new Parameter[childNum];
        for(int i = 0; i < modelIndicators.length; i++){
            modelIndicators[i] = (Parameter)xoInd.getChild(i);
        }
             return new MicrosatelliteModelSelectOperator(modelChoose, modelIndicators, weight);
    }
         //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************
         public String getParserDescription() {
        return "This element returns a microsatellite averaging operator on a given parameter.";
    }
         public Class getReturnType() {
        return MCMCOperator.class;
    }
         public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }
         private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
            new ElementRule(MODEL_CHOOSE, new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new ElementRule(MODEL_INDICATORS, new XMLSyntaxRule[]{new ElementRule(Parameter.class,1,Integer.MAX_VALUE)}),
    };
}