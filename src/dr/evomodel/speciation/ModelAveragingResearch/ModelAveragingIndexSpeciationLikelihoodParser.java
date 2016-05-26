/*
 * ModelAveragingIndexSpeciationLikelihoodParser.java
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

package dr.evomodel.speciation.ModelAveragingResearch;


import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.xml.*;


/**
 * @author Alexei Drummond
 * @author Walter Xie
 */
public class ModelAveragingIndexSpeciationLikelihoodParser extends AbstractXMLObjectParser {
    public static final String MODEL_AVE_Index_SPECIATION_LIKELIHOOD = "modelAveragingIndexSpeciationLikelihood";
    public static final String INDEX = "modelIndex";
    public static final String MAX_INDEX = "maxIndex";

    public String getParserName() {
        return MODEL_AVE_Index_SPECIATION_LIKELIHOOD;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Variable<Integer> index;

//        cxo = xo.getChild(INDEX);
        index = (Variable<Integer>) xo.getElementFirstChild(INDEX); // integer index parameter size = real size - 1
        Parameter maxIndex = (Parameter) xo.getElementFirstChild(MAX_INDEX);

        return new ModelAveragingIndexSpeciationLikelihood(xo.getId(), index, maxIndex);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "Model Averaging Speciation Likelihood.";
    }

    public Class getReturnType() {
        return ModelAveragingIndexSpeciationLikelihood.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(INDEX, new XMLSyntaxRule[]{
                    new ElementRule(Variable.class)
            }),
            new ElementRule(MAX_INDEX, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }, true),
    };

}
