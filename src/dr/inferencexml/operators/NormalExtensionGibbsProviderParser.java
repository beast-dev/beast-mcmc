/*
 * NormalExtensionGibbsProviderParser.java
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

package dr.inferencexml.operators;

import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.preorder.ModelExtensionProvider;
import dr.inference.operators.repeatedMeasures.GammaGibbsProvider;
import dr.xml.*;

/**
 * @author Gabriel Hassler
 * @author Marc A. Suchard
 */

public class NormalExtensionGibbsProviderParser extends AbstractXMLObjectParser {

    private static final String NORMAL_EXTENSION = "normalExtension";
    private static final String TREE_TRAIT = "treeTraitName";


    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        ModelExtensionProvider.NormalExtensionProvider dataModel = (ModelExtensionProvider.NormalExtensionProvider)
                xo.getChild(ModelExtensionProvider.NormalExtensionProvider.class);

        TreeDataLikelihood likelihood = (TreeDataLikelihood) xo.getChild(TreeDataLikelihood.class);

        String traitName = null;
        if (xo.hasAttribute(TREE_TRAIT)) {
            traitName = xo.getStringAttribute(TREE_TRAIT);
        }

        return new GammaGibbsProvider.NormalExtensionGibbsProvider(dataModel, likelihood, traitName);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                new ElementRule(ModelExtensionProvider.NormalExtensionProvider.class),
                new ElementRule(TreeDataLikelihood.class),
                AttributeRule.newStringRule(TREE_TRAIT, true)
        };
    }

    @Override
    public String getParserDescription() {
        return "Provides sufficient statistics for normal precision with gamma prior for extended tree model.";
    }

    @Override
    public Class getReturnType() {
        return GammaGibbsProvider.NormalExtensionGibbsProvider.class;
    }

    @Override
    public String getParserName() {
        return NORMAL_EXTENSION;
    }
}
