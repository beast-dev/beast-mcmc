/*
 * LatentFactorModelParser.java
 *
 * Copyright (c) 2002-2013 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

import dr.evomodel.continuous.LatentFactorModel;
import dr.evomodel.tree.TreeModel;
import dr.evomodelxml.treelikelihood.TreeTraitParserUtilities;
import dr.inference.model.CompoundParameter;
import dr.inference.model.MatrixParameter;
import dr.xml.*;

import java.util.List;

/**
 * @author Max Tolkoff
 * @author Marc Suchard
 */

public class LatentFactorModelParser extends AbstractXMLObjectParser {
    public final static String LATENT_FACTOR_MODEL = "latentFactorModel";
    public final static String NUMBER_OF_FACTORS = "factorNumber";
    public final static String FACTORS = "factors";
    public final static String DATA = "data";
    public final static String LOADINGS = "loadings";
    public static final String PRECISION = "precision";

    public String getParserName() {
        return LATENT_FACTOR_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
//        Parameter latent  = null;
        CompoundParameter factors = (CompoundParameter) xo.getChild(FACTORS).getChild(CompoundParameter.class);

        TreeTraitParserUtilities utilities = new TreeTraitParserUtilities();
        String traitName = TreeTraitParserUtilities.DEFAULT_TRAIT_NAME;

        TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);

        TreeTraitParserUtilities.TraitsAndMissingIndices returnValue =
                utilities.parseTraitsFromTaxonAttributes(xo, traitName, treeModel, true);
        CompoundParameter dataParameter = returnValue.traitParameter;
        List<Integer> missingIndices = returnValue.missingIndices;
        traitName = returnValue.traitName;
//
//
//
//
//        MatrixParameter data = (MatrixParameter) xo.getChild(DATA).getChild(MatrixParameter.class);
        MatrixParameter loadings = (MatrixParameter) xo.getChild(LOADINGS).getChild(MatrixParameter.class);
        MatrixParameter precision = (MatrixParameter) xo.getChild(PRECISION).getChild(MatrixParameter.class);
        int numFactors = xo.getAttribute(NUMBER_OF_FACTORS, 4);


        return new LatentFactorModel(dataParameter, factors, loadings, precision, numFactors);
    }

    private static final XMLSyntaxRule[] rules = {
            AttributeRule.newIntegerRule(NUMBER_OF_FACTORS),
            new ElementRule(TreeModel.class),
            new ElementRule(FACTORS, new XMLSyntaxRule[]{
                    new ElementRule(CompoundParameter.class),
            }),
            new ElementRule(DATA, new XMLSyntaxRule[]{
                    new ElementRule(CompoundParameter.class)
            }),
            new ElementRule(LOADINGS, new XMLSyntaxRule[]{
                    new ElementRule(MatrixParameter.class)
            }),
            new ElementRule(PRECISION, new XMLSyntaxRule[]{
                    new ElementRule(MatrixParameter.class)
            }),
    };

//    <latentFactorModel>
//      <factors>
//         <parameter idref="factors"/>
//      </factors>
//    </latentFactorModel>


    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    @Override
    public String getParserDescription() {
        return "Sets up a latent factor model, with starting guesses for the loadings and factor matrices as well as the data for the factor analysis";
    }

    @Override
    public Class getReturnType() {
        return LatentFactorModel.class;
    }


}
