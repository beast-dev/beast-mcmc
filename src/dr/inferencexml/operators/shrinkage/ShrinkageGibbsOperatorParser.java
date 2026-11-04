/*
 * ShrinkageGibbsOperatorParser.java
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

import dr.inference.distribution.IndependentInverseGammaDistributionModel;
import dr.inference.operators.shrinkage.ShrinkageGibbsOperator;
import dr.inference.model.Parameter;
import dr.xml.*;

public class ShrinkageGibbsOperatorParser extends AbstractXMLObjectParser{
    private final static String SHRINKAGE_GIBBS_OPERATOR = "shrinkageGibbsOperator";
    private final static String LOCAL_PRIOR = "localPrior";
    private final static String GLOBAL_PRIOR = "globalPrior";
    private final static String DATA = "data";
    private final static String WEIGHT = "weight";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        IndependentInverseGammaDistributionModel localPrior = (IndependentInverseGammaDistributionModel) xo.getChild(LOCAL_PRIOR).getChild(IndependentInverseGammaDistributionModel.class);
        IndependentInverseGammaDistributionModel globalPrior = (IndependentInverseGammaDistributionModel) xo.getChild(GLOBAL_PRIOR).getChild(IndependentInverseGammaDistributionModel.class);
        Parameter data = (Parameter) xo.getChild(DATA).getChild(Parameter.class);

        double weight = xo.getDoubleAttribute(WEIGHT);


        return new ShrinkageGibbsOperator(weight, localPrior, globalPrior, data);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            AttributeRule.newDoubleRule(WEIGHT),
            new ElementRule(LOCAL_PRIOR, new XMLSyntaxRule[]{
                    new ElementRule(IndependentInverseGammaDistributionModel.class)
            }),
            new ElementRule(GLOBAL_PRIOR, new XMLSyntaxRule[]{
                    new ElementRule(IndependentInverseGammaDistributionModel.class)
            }),
            new ElementRule(DATA, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            })
    };

    @Override
    public String getParserDescription() {
        return null;
    }

    @Override
    public Class getReturnType() {
        return ShrinkageGibbsOperator.class;
    }

    @Override
    public String getParserName() {
        return SHRINKAGE_GIBBS_OPERATOR;
    }
}
