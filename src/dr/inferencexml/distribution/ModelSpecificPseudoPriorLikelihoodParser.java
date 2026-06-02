/*
 * ModelSpecificPseudoPriorLikelihoodParser.java
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

package dr.inferencexml.distribution;

import dr.xml.AbstractXMLObjectParser;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;
import dr.inference.distribution.DistributionLikelihood;
import dr.inference.distribution.ModelSpecificPseudoPriorLikelihood;
import dr.inference.model.Parameter;
import dr.math.distributions.Distribution;

/**
 * @author Chieh-Hsi Wu
 */
public class ModelSpecificPseudoPriorLikelihoodParser  extends AbstractXMLObjectParser {
    public static final String MODEL_SPECIFIC_PSEUDO_PRIOR = "modelSpecificPseudoPrior";
    public static final String PRIOR = "priorLik";
    public static final String PSEUDO_PRIOR = "pseudoPriorLik";
    public static final String MODELS = "models";
    public static final String MODEL_INDICATOR = "modelIndicator";
    public static final String SELECTED_VARIABLE = "selectedVariable";


    public String getParserName() {
        return MODEL_SPECIFIC_PSEUDO_PRIOR;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        DistributionLikelihood priorLikelihood = (DistributionLikelihood)xo.getElementFirstChild(PRIOR);
        DistributionLikelihood pseudoPriorLikelihood = (DistributionLikelihood)xo.getElementFirstChild(PSEUDO_PRIOR);
        Distribution prior = priorLikelihood.getDistribution();
        Distribution pseudoPrior = pseudoPriorLikelihood.getDistribution();
        Parameter modelIndicator = (Parameter)xo.getElementFirstChild(MODEL_INDICATOR);
        int[] models = xo.getIntegerArrayAttribute(MODELS);
        Parameter selectedVariable = (Parameter)xo.getElementFirstChild(SELECTED_VARIABLE);


        ModelSpecificPseudoPriorLikelihood likelihood =
                new ModelSpecificPseudoPriorLikelihood(
                        prior,
                        pseudoPrior,
                        modelIndicator,
                        models
                );
        likelihood.addData(selectedVariable);
        

        return likelihood;
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {

    };

    public String getParserDescription() {
        return "Calculates the likelihood of some data given some parametric or empirical distribution.";
    }

    public Class getReturnType() {
        return ModelSpecificPseudoPriorLikelihood.class;
    }
}
