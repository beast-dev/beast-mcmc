/*
 * ContinuousTraitDataModelParser.java
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

package dr.evomodelxml.continuous;

import dr.evolution.tree.Tree;
import dr.evomodel.treedatalikelihood.continuous.ContinuousTraitDataModel;
import dr.evomodel.treedatalikelihood.continuous.HomogeneousDiffusionModelDelegate;
import dr.evomodel.treedatalikelihood.continuous.IntegratedProcessTraitDataModel;
import dr.evomodel.treedatalikelihood.continuous.cdi.PrecisionType;
import dr.evomodelxml.treelikelihood.TreeTraitParserUtilities;
import dr.inference.model.CompoundParameter;
import dr.inference.model.Parameter;
import dr.xml.*;


public class ContinuousTraitDataModelParser extends AbstractXMLObjectParser {

    private static String CONTINUOUS_TRAITS = "continuousTraitDataModel";

    public static final String INTEGRATED_PROCESS = "integratedProcess";
    public static final String FORCE_COMPLETELY_MISSING = "forceCompletelyMissing";
    public static final String FORCE_FULL_PRECISION = "forceFullPrecision";


    public static final String NUM_TRAITS = "numTraits";


    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        return parseContinuousTraitDataModel(xo);
    }

    public static ContinuousTraitDataModel parseContinuousTraitDataModel(XMLObject xo) throws XMLParseException {
        return parseContinuousTraitDataModel(xo, null);
    }

    public static ContinuousTraitDataModel parseContinuousTraitDataModel(XMLObject xo, PrecisionType precisionType) throws XMLParseException {
        Tree treeModel = (Tree) xo.getChild(Tree.class);
        boolean[] missingIndicators;
        final String traitName;

        boolean useMissingIndices;
        boolean integratedProcess = xo.getAttribute(INTEGRATED_PROCESS, false);


        TreeTraitParserUtilities utilities = new TreeTraitParserUtilities();

        TreeTraitParserUtilities.TraitsAndMissingIndices returnValue =
                utilities.parseTraitsFromTaxonAttributes(xo, treeModel, true);
        CompoundParameter traitParameter = returnValue.traitParameter;

        int dimAll = traitParameter.getParameter(0).getDimension();
        int numTraits = xo.getAttribute(NUM_TRAITS, 1);
        int dim = dimAll / numTraits; //TODO: check that dimAll is a factor of numTraits, also TODO: maybe pass numTraits directly?


        missingIndicators = returnValue.getMissingIndicators();
//            sampleMissingParameter = returnValue.sampleMissingParameter;
        traitName = returnValue.traitName;
        useMissingIndices = returnValue.useMissingIndices;

        if (precisionType == null) {
            precisionType = PrecisionType.SCALAR;

            if (xo.getAttribute(FORCE_FULL_PRECISION, false) ||
                    (useMissingIndices && !xo.getAttribute(FORCE_COMPLETELY_MISSING, false))
//                    || !(diffusionProcessDelegate instanceof HomogeneousDiffusionModelDelegate)
                    ) {
                precisionType = PrecisionType.FULL;
            }
        }


        if (xo.hasChildNamed(TreeTraitParserUtilities.JITTER)) {
            utilities.jitter(xo, dim, missingIndicators);
        }

//            System.err.println("Using precisionType == " + precisionType + " for data model.");

        if (integratedProcess) {
            return new IntegratedProcessTraitDataModel(traitName,
                    traitParameter,
                    missingIndicators, useMissingIndices,
                    dim, precisionType);
        }

        return new ContinuousTraitDataModel(traitName,
                traitParameter,
                missingIndicators, useMissingIndices,
                dim, precisionType);
    }

    public static final XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new ElementRule(Tree.class),
            new ElementRule(TreeTraitParserUtilities.TRAIT_PARAMETER, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }),
            AttributeRule.newBooleanRule(INTEGRATED_PROCESS, true),
            AttributeRule.newIntegerRule(NUM_TRAITS, true),
            AttributeRule.newBooleanRule(FORCE_COMPLETELY_MISSING, true),
            AttributeRule.newStringRule(TreeTraitParserUtilities.TRAIT_NAME, true)
    };

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    @Override
    public String getParserDescription() {
        return "parses continuous traits from a tree";
    }

    @Override
    public Class getReturnType() {
        return ContinuousTraitDataModel.class;
    }

    @Override
    public String getParserName() {
        return CONTINUOUS_TRAITS;
    }
}
