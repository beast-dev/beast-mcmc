/*
 * VarianceProportionStatisticParser.java
 *
 * Copyright (c) 2002-2020 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.inferencexml.model;

import dr.evolution.tree.Tree;
import dr.evomodel.continuous.MultivariateDiffusionModel;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.continuous.ContinuousDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.continuous.ContinuousTraitPartialsProvider;
import dr.evomodel.treedatalikelihood.continuous.RepeatedMeasuresTraitDataModel;
import dr.evomodel.treedatalikelihood.continuous.TreeScaledRepeatedMeasuresTraitDataModel;
import dr.evomodelxml.treedatalikelihood.ContinuousDataLikelihoodParser;
import dr.inference.model.VarianceProportionStatistic;
import dr.inference.model.VarianceProportionStatisticEmpirical;
import dr.inference.model.VarianceProportionStatisticPopulation;
import dr.xml.*;

/**
 * @author Gabriel Hassler
 */

public class VarianceProportionStatisticParser extends AbstractXMLObjectParser {

    public static final String PARSER_NAME = "varianceProportionStatistic";
    private static final String MATRIX_RATIO = "matrixRatio";
    private static final String ELEMENTWISE = "elementWise";
    private static final String SYMMETRIC_DIVISION = "symmetricDivision";
    private static final String CO_HERITABILITY = "coheritability";
    private static final String EMPIRICAL = "useEmpiricalVariance";
    private static final String POPULATION = "usePopulationVariance";
    private static final String FORCE_SAMPLING = "forceSampling";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        TreeDataLikelihood treeLikelihood = (TreeDataLikelihood) xo.getChild(TreeDataLikelihood.class);

        Tree tree = treeLikelihood.getTree();

        ContinuousTraitPartialsProvider dataModel = ((ContinuousDataLikelihoodDelegate) treeLikelihood.getDataLikelihoodDelegate()).getDataModel();
        if (!(dataModel instanceof RepeatedMeasuresTraitDataModel)) {
            throw new RuntimeException(
                    "In " + PARSER_NAME + ": " +
                            "The provided likelihood does not have a " + RepeatedMeasuresTraitDataModel.REPEATED_MEASURES_MODEL + " element." +
                            "VarianceProportionStatistic is only implemented for repeated measures.");
        }

        boolean empirical = xo.getAttribute(EMPIRICAL, false);
        boolean forceSampling = xo.getAttribute(FORCE_SAMPLING, true);
        boolean population = xo.getAttribute(POPULATION, false);

        if (dataModel instanceof TreeScaledRepeatedMeasuresTraitDataModel && !population) {
            throw new RuntimeException(
                    "varianceProportionStatistic with " +
                            POPULATION + "=false" + " is not yet implemented for " +
                            "repeatedMeasuresModel argument scaleByTreeHeight='true'.");
        }

        MultivariateDiffusionModel diffusionModel = ((ContinuousDataLikelihoodDelegate) treeLikelihood.getDataLikelihoodDelegate()).getDiffusionModel();

        // If provided, check that the tree, data model and diffusion are consistent (backward compatibility)
        Tree treeXML = (Tree) xo.getChild(TreeModel.class);
        if ((treeXML != null) && (tree != treeXML)) {
            throw new RuntimeException(
                    "In " + PARSER_NAME + ": " +
                            "The provided tree is different from the tree in object " + ContinuousDataLikelihoodParser.CONTINUOUS_DATA_LIKELIHOOD + ".");
        }
        RepeatedMeasuresTraitDataModel dataModelXMl = (RepeatedMeasuresTraitDataModel)
                xo.getChild(RepeatedMeasuresTraitDataModel.class);
        if ((dataModelXMl != null) && (dataModel != dataModelXMl)) {
            throw new RuntimeException(
                    "In " + PARSER_NAME + ": " +
                            "The provided data model is different from the data model in object " + ContinuousDataLikelihoodParser.CONTINUOUS_DATA_LIKELIHOOD + ".");
        }
        MultivariateDiffusionModel diffusionModelXML = (MultivariateDiffusionModel)
                xo.getChild(MultivariateDiffusionModel.class);
        if ((diffusionModelXML != null) && (diffusionModel != diffusionModelXML)) {
            throw new RuntimeException(
                    "In " + PARSER_NAME + ": " +
                            "The provided diffusion model is different from the diffusion model in object " + ContinuousDataLikelihoodParser.CONTINUOUS_DATA_LIKELIHOOD + ".");
        }

        String ratioString = xo.getStringAttribute(MATRIX_RATIO);

        VarianceProportionStatistic.MatrixRatios ratio = null;

        if (ratioString.equalsIgnoreCase(ELEMENTWISE)) {
            ratio = VarianceProportionStatistic.MatrixRatios.ELEMENT_WISE;
        } else if (ratioString.equalsIgnoreCase(SYMMETRIC_DIVISION)) {
            ratio = VarianceProportionStatistic.MatrixRatios.SYMMETRIC_DIVISION;
        } else if (ratioString.equalsIgnoreCase(CO_HERITABILITY)) {
            ratio = VarianceProportionStatistic.MatrixRatios.CO_HERITABILITY;
        } else {
            throw new RuntimeException(PARSER_NAME + " must have attibute " + MATRIX_RATIO +
                    " with one of the following values: " + VarianceProportionStatistic.MatrixRatios.values());
        }

        if (empirical && population) {
            throw new RuntimeException(PARSER_NAME + "cannot use both empirical and population variances. Please set one to false.");
        }

        if (empirical) {
            return new VarianceProportionStatisticEmpirical(tree, treeLikelihood,
                    (RepeatedMeasuresTraitDataModel) dataModel, diffusionModel,
                    ratio, forceSampling);
        } else if (population) {
            return new VarianceProportionStatisticPopulation(tree, treeLikelihood,
                    (RepeatedMeasuresTraitDataModel) dataModel, diffusionModel, ratio);
        }
        return new VarianceProportionStatistic(tree, treeLikelihood,
                (RepeatedMeasuresTraitDataModel) dataModel, diffusionModel, ratio);
    }

    private final XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            AttributeRule.newStringRule(MATRIX_RATIO, false),
            AttributeRule.newStringRule(FORCE_SAMPLING, true),
            AttributeRule.newStringRule(EMPIRICAL, true),
            new ElementRule(TreeModel.class, true),
            new ElementRule(TreeDataLikelihood.class),
            new ElementRule(RepeatedMeasuresTraitDataModel.class, true),
            new ElementRule(MultivariateDiffusionModel.class, true)
    };

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    @Override
    public String getParserDescription() {
        return "This element returns a statistic that computes proportion of variance due to diffusion on the tree";
    }

    @Override
    public Class getReturnType() {
        return VarianceProportionStatistic.class;
    }

    @Override
    public String getParserName() {
        return PARSER_NAME;
    }
}
