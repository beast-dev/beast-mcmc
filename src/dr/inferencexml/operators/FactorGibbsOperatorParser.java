/*
 * FactorGibbsOperatorParser.java
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

import dr.inference.model.DiagonalMatrix;
import dr.inference.model.LatentFactorModel;
import dr.inference.operators.FactorGibbsOperator;
import dr.xml.*;

/**
 * @author Max Tolkoff
 */
public class FactorGibbsOperatorParser extends AbstractXMLObjectParser {
    private final String FACTOR_GIBBS_SAMPLER = "factorGibbsOperator";
    private final String WEIGHT = "weight";
    private final String RANDOM_SCAN = "randomScan";


    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        String weightTemp = (String) xo.getAttribute(WEIGHT);
        double weight = Double.parseDouble(weightTemp);
        DiagonalMatrix diffusionMatrix;
        diffusionMatrix = (DiagonalMatrix) xo.getChild(DiagonalMatrix.class);
        LatentFactorModel LFM = (LatentFactorModel) xo.getChild(LatentFactorModel.class);
        boolean randomScan = xo.getAttribute(RANDOM_SCAN, true);
        return new FactorGibbsOperator(LFM, weight, randomScan, diffusionMatrix);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new ElementRule(LatentFactorModel.class),
//            new ElementRule(CompoundParameter.class),
            new ElementRule(DiagonalMatrix.class),
            AttributeRule.newDoubleRule(WEIGHT),
    };

    @Override
    public String getParserDescription() {
        return null;
    }

    @Override
    public Class getReturnType() {
        return FactorGibbsOperator.class;
    }

    @Override
    public String getParserName() {
        return FACTOR_GIBBS_SAMPLER;
    }
}
