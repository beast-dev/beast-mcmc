/*
 * DataFromTreeTipsParser.java
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

package dr.evomodelxml.continuous;

import dr.evolution.tree.MultivariateTraitTree;
import dr.evomodelxml.treelikelihood.TreeTraitParserUtilities;
import dr.inference.model.MatrixParameter;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * @author Max Tolkoff
 * @author Marc Suchard
 */

public class DataFromTreeTipsParser extends AbstractXMLObjectParser {
    public final static String DATA_FROM_TREE_TIPS = "dataFromTreeTips";
    public final static String DATA = "data";
    public static final String CONTINUOUS = "continuous";
    public static final String SAMPLE_MISSING = "sampleMissing";
    public static final String MISSING_PARAMETER = "missingParameter";


    public String getParserName() {
        return DATA_FROM_TREE_TIPS;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        TreeTraitParserUtilities utilities = new TreeTraitParserUtilities();
        String traitName = (String) xo.getAttribute(TreeTraitParserUtilities.TRAIT_NAME);

        MultivariateTraitTree treeModel = (MultivariateTraitTree) xo.getChild(MultivariateTraitTree.class);
        boolean sampleMissing;
        if(xo.hasAttribute(SAMPLE_MISSING))
            sampleMissing = xo.getBooleanAttribute(SAMPLE_MISSING);
        else
            sampleMissing = true;

        TreeTraitParserUtilities.TraitsAndMissingIndices returnValue =
                utilities.parseTraitsFromTaxonAttributes(xo, traitName, treeModel, sampleMissing);
        MatrixParameter dataParameter = MatrixParameter.recast(returnValue.traitParameter.getId(),
                returnValue.traitParameter);

        if (xo.hasChildNamed(MISSING_PARAMETER)) {
            Parameter missing = (Parameter) xo.getChild(MISSING_PARAMETER).getChild(Parameter.class);
            missing.setDimension(dataParameter.getDimension());

            for (int i = 0; i < missing.getDimension(); i++) {
                if (returnValue.missingIndices.contains(i)) {
                    missing.setParameterValue(i, 1);
                } else {
                    missing.setParameterValue(i, 0);
                }
            }
        }

        return dataParameter;
    }

    private static final XMLSyntaxRule[] rules = {
            new ElementRule(MultivariateTraitTree.class),
            AttributeRule.newStringRule(TreeTraitParserUtilities.TRAIT_NAME),
            AttributeRule.newBooleanRule(SAMPLE_MISSING, true),
            new ElementRule(TreeTraitParserUtilities.TRAIT_PARAMETER, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }),
            new ElementRule(TreeTraitParserUtilities.MISSING, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }, true),
    };

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    @Override
    public String getParserDescription() {
        return "Takes the data from the tips of a tree and puts it into a MatrixParameter";
    }

    @Override
    public Class getReturnType() {
        return MatrixParameter.class;
    }


}
