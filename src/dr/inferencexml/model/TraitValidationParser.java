/*
 * TraitValidationParser.java
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

package dr.inferencexml.model;

import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodelxml.treelikelihood.TreeTraitParserUtilities;
import dr.inference.model.CrossValidationProvider;
import dr.inference.model.Parameter;
import dr.inference.model.TraitValidationProvider;
import dr.xml.*;

/**
 * @author Gabriel Hassler
 */

@Deprecated
public class TraitValidationParser extends AbstractXMLObjectParser {

    private static final String TRAIT_VALIDATION = "traitValidation";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        TraitValidationProvider provider = TraitValidationProviderParser.parseTraitValidationProvider(xo);

        boolean logSum = xo.getAttribute(CrossValidatorParser.LOG_SUM, false);

        if (logSum)
            return new CrossValidationProvider.CrossValidatorSum(provider, CrossValidationProvider.ValidationType.SQUARED_ERROR);
        return new CrossValidationProvider.CrossValidator(provider, CrossValidationProvider.ValidationType.SQUARED_ERROR);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                new ElementRule(TreeDataLikelihood.class),
                AttributeRule.newStringRule(TreeTraitParserUtilities.TRAIT_NAME),
                AttributeRule.newStringRule(TraitValidationProviderParser.INFERRED_NAME),
                new ElementRule(TreeTraitParserUtilities.TRAIT_PARAMETER, new XMLSyntaxRule[]{
                        new ElementRule(Parameter.class)
                }),
                new ElementRule(TraitValidationProviderParser.MASK, new XMLSyntaxRule[]{
                        new ElementRule(Parameter.class)
                }, true),
                AttributeRule.newBooleanRule(CrossValidatorParser.LOG_SUM)
        };
    }

    @Override
    public String getParserDescription() {
        return null;
    }

    @Override
    public Class getReturnType() {
        return CrossValidationProvider.CrossValidator.class;
    }

    @Override
    public String getParserName() {
        return TRAIT_VALIDATION;
    }
}
