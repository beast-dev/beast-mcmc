/*
 * ShrinkageAugmentedGibbsOperatorParser.java
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

package dr.inferencexml.operators.shrinkage;

import dr.inference.distribution.DistributionLikelihood;
import dr.inference.distribution.IndependentInverseGammaDistributionModel;
import dr.inference.operators.shrinkage.ShrinkageGibbsOperator;
import dr.inference.operators.shrinkage.ShrinkageAugmentedGibbsOperator;
import dr.xml.*;

public class ShrinkageAugmentedGibbsOperatorParser extends AbstractXMLObjectParser{
    public final static String LOCAL = "local";
    public final static String GLOBAL = "global";
    public final static String SHRINKAGE_AUGMENTED_GIBBS_OPERATOR = "shrinkageAugmentedGibbsOperator";
    public final static String WEIGHT = "weight";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        double weight = xo.getDoubleAttribute(WEIGHT);
        XMLObject local = xo.getChild(LOCAL);
        XMLObject global = xo.getChild(GLOBAL);

        IndependentInverseGammaDistributionModel localPrior = (IndependentInverseGammaDistributionModel) local.getChild(IndependentInverseGammaDistributionModel.class);
        DistributionLikelihood localAugmentedPrior = (DistributionLikelihood) local.getChild(DistributionLikelihood.class);

        IndependentInverseGammaDistributionModel globalPrior = (IndependentInverseGammaDistributionModel) global.getChild(IndependentInverseGammaDistributionModel.class);
        DistributionLikelihood globalAugmentedPrior = (DistributionLikelihood) global.getChild(DistributionLikelihood.class);



        return new ShrinkageAugmentedGibbsOperator(weight, localAugmentedPrior, globalAugmentedPrior, localPrior, globalPrior);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            AttributeRule.newDoubleRule(WEIGHT),
            new ElementRule(LOCAL, new XMLSyntaxRule[]{
                    new ElementRule(IndependentInverseGammaDistributionModel.class),
                    new ElementRule(DistributionLikelihood.class)
            }),
            new ElementRule(GLOBAL, new XMLSyntaxRule[]{
                    new ElementRule(IndependentInverseGammaDistributionModel.class),
                    new ElementRule(DistributionLikelihood.class)
            }),
    };

    @Override
    public String getParserDescription() {
        return "Gibbs sampler for augmented portions of a shrinkage prior";
    }

    @Override
    public Class getReturnType() {
        return ShrinkageGibbsOperator.class;
    }

    @Override
    public String getParserName() {
        return SHRINKAGE_AUGMENTED_GIBBS_OPERATOR;
    }
}
