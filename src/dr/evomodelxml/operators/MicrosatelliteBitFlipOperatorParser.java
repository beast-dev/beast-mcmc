/*
 * MsatBitFlipOperatorParser.java
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

import dr.xml.*;
import dr.inference.model.Parameter;
import dr.inference.operators.MCMCOperator;
import dr.evomodel.operators.MicrosatelliteBitFlipOperator;

/**
 * @author Chieh-Hsi Wu
 *
 * Parser for MicrosatelliteBitFlipOperator
 */
public class MicrosatelliteBitFlipOperatorParser extends AbstractXMLObjectParser{
    public static final String MODEL_CHOOSE = "modelChoose";
    public static final String DEPENDENCIES = "dependencies";
    public static final String VARIABLE_INDICES = "variableIndices";

    public String getParserName() {
        return "msatModelSwitchOperator";
    }
         public Object parseXMLObject(XMLObject xo) throws XMLParseException {
             double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);
        Parameter modelChoose = (Parameter) xo.getElementFirstChild(MODEL_CHOOSE);
        Parameter dependencies = (Parameter)xo.getElementFirstChild(DEPENDENCIES);
        int[] variableIndices;
            if(xo.hasChildNamed(VARIABLE_INDICES)){

                double[] temp = ((Parameter)xo.getElementFirstChild(VARIABLE_INDICES)).getParameterValues();
                variableIndices = new int[temp.length];
                for(int i = 0; i < temp.length;i++){
                    variableIndices[i] = (int)temp[i];
                }

            }else{
                variableIndices = new int[]{0, 1, 2, 3, 4, 5};
            }

            return new MicrosatelliteBitFlipOperator(modelChoose, dependencies, weight, variableIndices);
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
            new ElementRule(DEPENDENCIES, new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new ElementRule(VARIABLE_INDICES, new XMLSyntaxRule[]{new ElementRule(Parameter.class)},true)
    };
}
