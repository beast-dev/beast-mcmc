/*
 * TwoPhaseModelParser.java
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
import dr.oldevomodel.substmodel.OnePhaseModel;
import dr.oldevomodel.substmodel.TwoPhaseModel;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * @author Chieh-Hsi Wu
 *
 * Parser for Two Phase Model of microsatellites.
 *
 */
public class TwoPhaseModelParser extends AbstractXMLObjectParser{
    public static final String SUBMODEL = "Submodel";
    public static final String GEO_PARAM = "GeoParam";
    public static final String ONEPHASEPR_PARAM ="OnePhasePrParam";
    public static final String TRANS_PARAM = "TransformParam";
    public static final double UPPER = 1.0;
    public static final double LOWER = 0.0;
    public static final String ESTIMATE_SUBMODEL_PARAMS = "estimateSubmodelParameters";

    public String getParserName() {
        return TwoPhaseModel.TWO_PHASE_MODEL;
    }

    //AbstractXMLObjectParser implementation
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        OnePhaseModel subModel = (OnePhaseModel) xo.getElementFirstChild(SUBMODEL);
        Microsatellite dataType = (Microsatellite)xo.getChild(Microsatellite.class);

        Parameter.Default geoParam =(Parameter.Default) xo.getElementFirstChild(GEO_PARAM);

        Parameter paramP = (Parameter) xo.getElementFirstChild(ONEPHASEPR_PARAM);

        Parameter limitE = null;
        if(xo.hasChildNamed(TRANS_PARAM)){
            limitE = (Parameter) xo.getElementFirstChild(TRANS_PARAM);
        }

        boolean estimateSubmodelParams =  xo.getAttribute(ESTIMATE_SUBMODEL_PARAMS,false);

        FrequencyModel freqModel = null;
        if(xo.hasChildNamed(FrequencyModelParser.FREQUENCIES)){
            freqModel = (FrequencyModel)xo.getElementFirstChild(FrequencyModelParser.FREQUENCIES);
        }

        return new TwoPhaseModel(dataType, freqModel, subModel, paramP, geoParam, limitE, estimateSubmodelParams);
    }

    public String getParserDescription() {
        return "This element represents an instance of the Two Phase Model of microsatellite evolution.";
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
        new ElementRule(Microsatellite.class),
        new ElementRule(FrequencyModelParser.FREQUENCIES, new XMLSyntaxRule[]{new ElementRule(FrequencyModel.class)},true),
        new ElementRule(SUBMODEL,new XMLSyntaxRule[]{new ElementRule(OnePhaseModel.class)}),
        new ElementRule(ONEPHASEPR_PARAM,new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
        new ElementRule(GEO_PARAM,new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
        new ElementRule(TRANS_PARAM,new XMLSyntaxRule[]{new ElementRule(Parameter.class)},true),
        new StringAttributeRule(ESTIMATE_SUBMODEL_PARAMS,"whether or not to esitmate the parameters of the submodel",true)
    };

    public Class getReturnType() {
        return TwoPhaseModel.class;
    }


}
