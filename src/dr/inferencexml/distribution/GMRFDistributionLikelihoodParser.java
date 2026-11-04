/*
 * GMRFDistributionLikelihoodParser.java
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

import dr.evomodelxml.coalescent.GMRFSkyrideLikelihoodParser;
import dr.evomodelxml.coalescent.smooth.SmoothSkygridLikelihoodParser;
import dr.inference.model.Parameter;
import dr.math.distributions.GmrfDistributionLikelihood;
import dr.xml.*;

/**
 * @author Marc Suchard
 * @author Xiang Ji
 */
public class GMRFDistributionLikelihoodParser  extends AbstractXMLObjectParser {

    private static final String DATA = "data";
    private static final String PRECISION_PARAMETER = GMRFSkyrideLikelihoodParser.PRECISION_PARAMETER;

    private static final String SINGLE_BETA = GMRFSkyrideLikelihoodParser.SINGLE_BETA;
    private static final String TIME_AWARE_SMOOTHING = GMRFSkyrideLikelihoodParser.TIME_AWARE_SMOOTHING;

    private static final String LAMBDA_PARAMETER = GMRFSkyrideLikelihoodParser.LAMBDA_PARAMETER;


    public final String NAME = "gmrfDistributionLikelihood";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        XMLObject cxo = xo.getChild(DATA);
        Parameter data = (Parameter) cxo.getChild(Parameter.class);

        cxo = xo.getChild(PRECISION_PARAMETER);
        Parameter precParameter = (Parameter) cxo.getChild(Parameter.class);

        Parameter lambda;
        if (xo.getChild(LAMBDA_PARAMETER) != null) {
            cxo = xo.getChild(LAMBDA_PARAMETER);
            lambda = (Parameter) cxo.getChild(Parameter.class);
        } else {
            lambda = new Parameter.Default(LAMBDA_PARAMETER, 1.0);
        }

        Parameter gridPoints = SmoothSkygridLikelihoodParser.getGridPoints(xo);
        return new GmrfDistributionLikelihood(NAME, precParameter, lambda, gridPoints, data);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }
    private final XMLSyntaxRule[] rules = {
            new ElementRule(DATA, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }),
            new ElementRule(PRECISION_PARAMETER, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }),
            new ElementRule(SINGLE_BETA, new XMLSyntaxRule[] {
                    new ElementRule(Parameter.class),
            }, true),
            AttributeRule.newBooleanRule(TIME_AWARE_SMOOTHING, true),
    };

    @Override
    public String getParserDescription() {
        return null;
    }

    @Override
    public Class getReturnType() {
        return GmrfDistributionLikelihood.class;
    }

    @Override
    public String getParserName() {
        return NAME;
    }
}
