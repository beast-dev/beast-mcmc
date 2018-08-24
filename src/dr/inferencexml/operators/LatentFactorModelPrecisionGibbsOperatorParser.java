/*
 * LatentFactorModelPrecisionGibbsOperatorParser.java
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

package dr.inferencexml.operators;

import dr.inference.distribution.DistributionLikelihood;
import dr.inference.model.LatentFactorModel;
import dr.inference.operators.LatentFactorModelPrecisionGibbsOperator;
import dr.xml.*;

/**
 * Created by max on 6/12/14.
 */
public class LatentFactorModelPrecisionGibbsOperatorParser extends AbstractXMLObjectParser {
    public final String LATENT_FACTOR_MODEL_PRECISION_OPERATOR = "latentFactorModelPrecisionOperator";
    public final String WEIGHT = "weight";
    public final String RANDOM_SCAN = "randomScan";
    public final String PATH_PARAMETER = "pathParameter";


    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        String weightTemp = (String) xo.getAttribute(WEIGHT);
        double weight = Double.parseDouble(weightTemp);
        LatentFactorModel LFM = (LatentFactorModel) xo.getChild(LatentFactorModel.class);
        DistributionLikelihood prior = (DistributionLikelihood) xo.getChild(DistributionLikelihood.class);
        boolean randomScan = xo.getAttribute(RANDOM_SCAN, true);
        LatentFactorModelPrecisionGibbsOperator lfmOp = new LatentFactorModelPrecisionGibbsOperator(LFM, prior, weight, randomScan);
        if(xo.hasAttribute(PATH_PARAMETER)){
            System.out.println("WARNING: Setting Path Parameter is intended for debugging purposes only!");
            lfmOp.setPathParameter(xo.getDoubleAttribute(PATH_PARAMETER));
        }

        return lfmOp;


    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new ElementRule(LatentFactorModel.class),
            new ElementRule(DistributionLikelihood.class),
//            new ElementRule(CompoundParameter.class),
            AttributeRule.newDoubleRule(WEIGHT),
            AttributeRule.newDoubleRule(PATH_PARAMETER, true),
    };

    @Override
    public String getParserDescription() {
        return "Gibbs sampler for the precision of a factor analysis model";
    }

    @Override
    public Class getReturnType() {
        return LatentFactorModelPrecisionGibbsOperator.class;
    }

    @Override
    public String getParserName() {
        return LATENT_FACTOR_MODEL_PRECISION_OPERATOR;
    }
}
