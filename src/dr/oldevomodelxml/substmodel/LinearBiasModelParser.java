/*
 * LinearBiasModelParser.java
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
import dr.oldevomodel.substmodel.FrequencyModel;
import dr.oldevomodel.substmodel.LinearBiasModel;
import dr.oldevomodel.substmodel.OnePhaseModel;
import dr.inference.model.Parameter;
import dr.xml.*;
/**
 * @author Chieh-Hsi Wu
 *
 * Parser for LinearBiasModel of Microsatellite.
 */
public class LinearBiasModelParser extends AbstractXMLObjectParser {
    public static final String SUBMODEL = "Submodel";
    public static final String BIAS_CONSTANT = "BiasConstant";
    public static final String BIAS_LINEAR = "BiasLinear";
    public static final String ESTIMATE_SUBMODEL_PARAMS = "estimateSubmodelParameters";
    public static final String LOGISTICS = "logistics";
    public static final String IS_SUBMODEL = "isSubmodel";


    public String getParserName() {
       return LinearBiasModel.LINEAR_BIAS_MODEL;
    }


    //AbstractXMLObjectParser implementation
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        OnePhaseModel subModel = (OnePhaseModel) xo.getElementFirstChild(SUBMODEL);
        Microsatellite dataType = (Microsatellite)subModel.getDataType();

        Parameter biasConst = null;
        if(xo.hasChildNamed(BIAS_CONSTANT)){
            biasConst =(Parameter) xo.getElementFirstChild(BIAS_CONSTANT);
        }

        Parameter biasLin = null;
        if(xo.hasChildNamed(BIAS_LINEAR)){
            biasLin = (Parameter) xo.getElementFirstChild(BIAS_LINEAR);
        }

        //get FrequencyModel
        FrequencyModel freqModel = null;
        if(xo.hasChildNamed(FrequencyModelParser.FREQUENCIES)){
            freqModel = (FrequencyModel)xo.getElementFirstChild(FrequencyModelParser.FREQUENCIES);
        }

        boolean estimateSubmodelParams = false;
        if(xo.hasAttribute(ESTIMATE_SUBMODEL_PARAMS)){
            estimateSubmodelParams = xo.getBooleanAttribute(ESTIMATE_SUBMODEL_PARAMS);
        }
        System.out.println("Is estimating submodel parameter(s): "+estimateSubmodelParams);

        boolean logistics = false;
        if(xo.hasAttribute(LOGISTICS)){
            logistics = xo.getBooleanAttribute(LOGISTICS);
        }
        System.out.println("Using logistic regression: "+ logistics);

        boolean isSubmodel = false;
        if(xo.hasAttribute(IS_SUBMODEL)){
            isSubmodel = xo.getBooleanAttribute(IS_SUBMODEL);
        }
        System.out.println("Is a submodel: "+isSubmodel);

        return new LinearBiasModel(
                dataType,
                freqModel,
                subModel,
                biasConst,
                biasLin,
                logistics,
                estimateSubmodelParams,
                isSubmodel
        );
    }

    public String getParserDescription() {
        return "This element represents an instance of the stepwise mutation model of microsatellite evolution.";
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new ElementRule(FrequencyModelParser.FREQUENCIES, new XMLSyntaxRule[]{
                    new ElementRule(FrequencyModel.class)},true),
            new ElementRule(SUBMODEL,new XMLSyntaxRule[]{new ElementRule(OnePhaseModel.class)}),
            new ElementRule(Microsatellite.class),
            new ElementRule(BIAS_CONSTANT,new XMLSyntaxRule[]{new ElementRule(Parameter.class)},true),
            new ElementRule(BIAS_LINEAR,new XMLSyntaxRule[]{new ElementRule(Parameter.class)},true),
            new StringAttributeRule(ESTIMATE_SUBMODEL_PARAMS,"whether or not to esitmate the parameters of the submodel",true),
            AttributeRule.newBooleanRule(LOGISTICS,true),
            AttributeRule.newBooleanRule(IS_SUBMODEL,true)
    };

    public Class getReturnType() {
        return LinearBiasModel.class;
    }


}
