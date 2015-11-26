/*
 * LoadingsIndependenceOperatorParser.java
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
import dr.inference.operators.CoercionMode;
import dr.inference.operators.LoadingsIndependenceOperator;
import dr.xml.*;

/**
 * Created with IntelliJ IDEA.
 * User: max
 * Date: 5/23/14
 * Time: 1:49 PM
 * To change this template use File | Settings | File Templates.
 */
public class LoadingsIndependenceOperatorParser extends AbstractXMLObjectParser {
    public static final String LOADINGS_GIBBS_OPERATOR = "loadingsIndependenceOperator";
    public static final String WEIGHT = "weight";
    private final String RANDOM_SCAN = "randomScan";
    private final String SCALE_FACTOR = "scaleFactor";


    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        CoercionMode mode = CoercionMode.parseMode(xo);
        String scaleFactorTemp = (String) xo.getAttribute(SCALE_FACTOR);
        double scaleFactor = Double.parseDouble(scaleFactorTemp);
        String weightTemp = (String) xo.getAttribute(WEIGHT);
        double weight = Double.parseDouble(weightTemp);
        LatentFactorModel LFM = (LatentFactorModel) xo.getChild(LatentFactorModel.class);
        DistributionLikelihood prior = (DistributionLikelihood) xo.getChild(DistributionLikelihood.class);
        boolean randomScan = xo.getAttribute(RANDOM_SCAN, true);

        return new LoadingsIndependenceOperator(LFM, prior, weight, randomScan, scaleFactor, mode);  //To change body of implemented methods use File | Settings | File Templates.
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
    };

    @Override
    public String getParserDescription() {
        return "Gibbs sampler for the loadings matrix of a latent factor model";  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public Class getReturnType() {
        return LoadingsIndependenceOperator.class;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public String getParserName() {
        return LOADINGS_GIBBS_OPERATOR;  //To change body of implemented methods use File | Settings | File Templates.
    }
}
