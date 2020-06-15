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
import dr.evomodel.treedatalikelihood.continuous.RepeatedMeasuresTraitDataModel;
import dr.evomodel.treedatalikelihood.continuous.TreeScaledRepeatedMeasuresTraitDataModel;
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

        Tree tree = (Tree) xo.getChild(TreeModel.class);
        RepeatedMeasuresTraitDataModel dataModel = (RepeatedMeasuresTraitDataModel)
                xo.getChild(RepeatedMeasuresTraitDataModel.class);

        boolean empirical = xo.getAttribute(EMPIRICAL, false);
        boolean forceSampling = xo.getAttribute(FORCE_SAMPLING, true);
        boolean population = xo.getAttribute(POPULATION, false);

        if (dataModel instanceof TreeScaledRepeatedMeasuresTraitDataModel && !population) {
            throw new RuntimeException(
                    "varianceProportionStatistic with " +
                            POPULATION + "=false" + " is not yet implemented for " +
                            "repeatedMeasuresModel argument scaleByTreeHeight='true'.");
        }

        MultivariateDiffusionModel diffusionModel = (MultivariateDiffusionModel)
                xo.getChild(MultivariateDiffusionModel.class);

        TreeDataLikelihood treeLikelihood = (TreeDataLikelihood) xo.getChild(TreeDataLikelihood.class);

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
            return new VarianceProportionStatisticEmpirical(tree, treeLikelihood, dataModel, diffusionModel,
                    ratio, forceSampling);
        } else if (population) {
            return new VarianceProportionStatisticPopulation(tree, treeLikelihood, dataModel, diffusionModel, ratio);
        }
        return new VarianceProportionStatistic(tree, treeLikelihood, dataModel, diffusionModel, ratio);
    }

    private final XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            AttributeRule.newStringRule(MATRIX_RATIO, false),
            AttributeRule.newStringRule(FORCE_SAMPLING, true),
            AttributeRule.newStringRule(EMPIRICAL, true),
            new ElementRule(TreeModel.class),
            new ElementRule(TreeDataLikelihood.class),
            new ElementRule(RepeatedMeasuresTraitDataModel.class),
            new ElementRule(MultivariateDiffusionModel.class)
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
