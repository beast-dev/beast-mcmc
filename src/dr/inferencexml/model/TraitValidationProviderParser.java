/*
 * TraitValidationProviderParser.java
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

import dr.evolution.tree.Tree;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.continuous.ContinuousDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.continuous.ContinuousTraitPartialsProvider;
import dr.evomodel.treedatalikelihood.preorder.ModelExtensionProvider;
import dr.evomodelxml.treelikelihood.TreeTraitParserUtilities;
import dr.inference.model.CrossValidationProvider;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.TraitValidationProvider;
import dr.xml.*;

/**
 * @author Gabriel Hassler
 */

public class TraitValidationProviderParser extends AbstractXMLObjectParser {

    public final static String TRAIT_VALIDATION_PROVIDER = "traitValidationProvider";
    final static String MASK = "mask";
    final static String INFERRED_NAME = "inferredTrait";


    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        return parseTraitValidationProvider(xo);

    }

    public static TraitValidationProvider parseTraitValidationProvider(XMLObject xo) throws XMLParseException {
        String inferredValuesName = xo.getStringAttribute(INFERRED_NAME);

        TreeDataLikelihood treeLikelihood = (TreeDataLikelihood) xo.getChild(TreeDataLikelihood.class);

        ContinuousDataLikelihoodDelegate delegate =
                (ContinuousDataLikelihoodDelegate) treeLikelihood.getDataLikelihoodDelegate();

        ContinuousTraitPartialsProvider dataModel = delegate.getDataModel();
        final ModelExtensionProvider extensionProvider;
        if (dataModel instanceof ModelExtensionProvider) {
            extensionProvider = (ModelExtensionProvider) dataModel;
        } else {
            throw new XMLParseException("Tree likelihood delegate does must implement ModelExtensionProvider");
        }

        Tree treeModel = treeLikelihood.getTree();


        TreeTraitParserUtilities utilities = new TreeTraitParserUtilities();


        TreeTraitParserUtilities.TraitsAndMissingIndices returnValue =
                utilities.parseTraitsFromTaxonAttributes(xo,
                        treeModel, true);

        Parameter trueParameter = returnValue.traitParameter;
        boolean[] trueMissing = returnValue.getMissingIndicators();
        Parameter missingParameter = null;
        if (xo.hasChildNamed(MASK)) {
            missingParameter = (Parameter) xo.getElementFirstChild(MASK);
        }


        String id = xo.getId();


        TraitValidationProvider provider = new TraitValidationProvider(trueParameter, extensionProvider, treeModel, id,
                missingParameter, treeLikelihood, inferredValuesName, trueMissing);

        return provider;
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{

                new ElementRule(TreeDataLikelihood.class),
                AttributeRule.newStringRule(TreeTraitParserUtilities.TRAIT_NAME),
                AttributeRule.newStringRule(INFERRED_NAME),
                new ElementRule(TreeTraitParserUtilities.TRAIT_PARAMETER, new XMLSyntaxRule[]{
                        new ElementRule(Parameter.class)
                }),
                new ElementRule(MASK, new XMLSyntaxRule[]{
                        new ElementRule(Parameter.class)
                }, true)
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
        return TRAIT_VALIDATION_PROVIDER;
    }


}