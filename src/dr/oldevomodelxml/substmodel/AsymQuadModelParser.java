/*
 * AsymQuadModelParser.java
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

import dr.evolution.datatype.Microsatellite;
import dr.oldevomodel.substmodel.AsymmetricQuadraticModel;
import dr.oldevomodel.substmodel.FrequencyModel;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * @author Chieh-Hsi Wu
 *
 * Parser for Asymmetric Quadratic Model
 */
public class AsymQuadModelParser extends AbstractXMLObjectParser{

    public static final String EXPANSION_CONSTANT = "ExpansionConstant";
    public static final String CONTRACTION_CONSTANT = "ContractionConstant";
    public static final String EXPANSION_LIN = "ExpansionLinear";
    public static final String CONTRACTION_LIN = "ContractionLinear";
    public static final String EXPANSION_QUAD = "ExpansionQuad";
    public static final String CONTRACTION_QUAD = "ContractionQuad";
    public static final String IS_SUBMODEL = "isSubmodel";

    public String getParserName() {
        return AsymmetricQuadraticModel.ASYMQUAD_MODEL;
    }


    //AbstractXMLObjectParser implementation
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Microsatellite microsatellite = (Microsatellite) xo.getChild(Microsatellite.class);

        Parameter expanConst = processModelParameter(xo, EXPANSION_CONSTANT);

        Parameter expanLin = processModelParameter(xo, EXPANSION_LIN);

        Parameter expanQuad = processModelParameter(xo, EXPANSION_QUAD);

        Parameter contractConst = processModelParameter(xo, CONTRACTION_CONSTANT);

        Parameter contractLin = processModelParameter(xo, CONTRACTION_LIN);

        Parameter contractQuad = processModelParameter(xo, CONTRACTION_QUAD);

        //get FrequencyModel
        FrequencyModel freqModel = null;
        if(xo.hasChildNamed(FrequencyModelParser.FREQUENCIES)){
            freqModel = (FrequencyModel)xo.getElementFirstChild(FrequencyModelParser.FREQUENCIES);
        }

        boolean isSubmodel = xo.getAttribute(IS_SUBMODEL, false);

        return new AsymmetricQuadraticModel(
                microsatellite,
                freqModel,
                expanConst,
                expanLin,
                expanQuad,
                contractConst,
                contractLin,
                contractQuad,
                isSubmodel
        );
    }

    private Parameter processModelParameter(XMLObject xo,
                                          String parameterName)throws XMLParseException{
        Parameter param = null;
        if(xo.hasChildNamed(parameterName)){
            XMLObject paramXO = xo.getChild(parameterName);
            param =(Parameter) paramXO.getChild(Parameter.class);           
        }
        return param;
    }

    public String getParserDescription() {
        return "This element represents an instance of the stepwise mutation model of microsatellite evolution.";
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new ElementRule(Microsatellite.class),
            new ElementRule(FrequencyModel.class,true),
            new ElementRule(EXPANSION_CONSTANT,new XMLSyntaxRule[]{new ElementRule(Parameter.class)},true),
            new ElementRule(CONTRACTION_CONSTANT,new XMLSyntaxRule[]{new ElementRule(Parameter.class)},true),
            new ElementRule(EXPANSION_LIN,new XMLSyntaxRule[]{new ElementRule(Parameter.class)},true),
            new ElementRule(CONTRACTION_LIN,new XMLSyntaxRule[]{new ElementRule(Parameter.class)},true),
            new ElementRule(EXPANSION_QUAD,new XMLSyntaxRule[]{new ElementRule(Parameter.class)},true),
            new ElementRule(CONTRACTION_QUAD,new XMLSyntaxRule[]{new ElementRule(Parameter.class)},true),
            AttributeRule.newBooleanRule(IS_SUBMODEL,true)
    };

    public Class getReturnType() {
        return AsymmetricQuadraticModel.class;
    }

    public static boolean requirePattern(){
        return true;
    }
}

