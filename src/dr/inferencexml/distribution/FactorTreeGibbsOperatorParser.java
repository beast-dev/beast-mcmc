/*
 * FactorTreeGibbsOperatorParser.java
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

import dr.evomodel.continuous.GibbsSampleFromTreeInterface;
import dr.inference.model.LatentFactorModel;
import dr.inference.operators.factorAnalysis.FactorTreeGibbsOperator;
import dr.xml.*;

/**
 * Created by max on 5/16/16.
 */
public class FactorTreeGibbsOperatorParser extends AbstractXMLObjectParser {
    public static final String FACTOR_TREE_GIBBS_OPERATOR_PARSER = "factorTreeGibbsOperator";
    public static final String WEIGHT = "weight";
    public static final String RANDOM_SCAN = "randomScan";
    public static final String WORKING_PRIOR = "workingPrior";
    public final String PATH_PARAMETER = "pathParameter";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        String weightTemp = (String) xo.getAttribute(WEIGHT);
        double weight = Double.parseDouble(weightTemp);
        LatentFactorModel lfm= (LatentFactorModel) xo.getChild(LatentFactorModel.class);
        GibbsSampleFromTreeInterface tree = (GibbsSampleFromTreeInterface) xo.getChild(GibbsSampleFromTreeInterface.class);
        GibbsSampleFromTreeInterface workingTree = null;
        if(xo.getChild(WORKING_PRIOR) != null){
            workingTree = (GibbsSampleFromTreeInterface) xo.getChild(WORKING_PRIOR).getChild(GibbsSampleFromTreeInterface.class);
        }
        boolean randomScan = xo.getAttribute(RANDOM_SCAN, true);

        FactorTreeGibbsOperator lfmOp = new FactorTreeGibbsOperator(weight, lfm, tree, randomScan);
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
            new ElementRule(GibbsSampleFromTreeInterface.class),
            AttributeRule.newDoubleRule(WEIGHT),
            AttributeRule.newBooleanRule(RANDOM_SCAN, true),
            new ElementRule(WORKING_PRIOR, new XMLSyntaxRule[]{
                    new ElementRule(GibbsSampleFromTreeInterface.class)
            }, true),
            AttributeRule.newDoubleRule(PATH_PARAMETER, true),
    };


    @Override
    public String getParserDescription() {
        return "Gibbs sample a factor row (tip) on a tree";
    }

    @Override
    public Class getReturnType() {
        return FactorTreeGibbsOperator.class;
    }

    @Override
    public String getParserName() {
        return FACTOR_TREE_GIBBS_OPERATOR_PARSER;
    }
}
