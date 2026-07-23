/*
 * BNPRSamplingLikelihoodParser.java
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

package dr.evomodelxml.coalescent;

import dr.evomodel.coalescent.BNPRSamplingLikelihood;
import dr.evomodel.coalescent.demographicmodel.DemographicModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.MatrixParameter;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * Created by mkarcher on 3/10/17.
 */
public class BNPRSamplingLikelihoodParser extends AbstractXMLObjectParser {
    public static final String SAMPLING_LIKELIHOOD = "bnprSamplingLikelihood";
    public static final String MODEL = "model";
    public static final String BETAS = "betas";
    public static final String POPULATION_TREE = "populationTree";
    public static final String EPOCH_WIDTHS = "epochWidths";
    public static final String WIDTHS = "widths";
    public static final String COVARIATES = "covariates";
    public static final String POWER_COVARIATES = "powerCovariates";
    public static final String POWER_BETAS = "powerBetas";

    @Override
    public String getParserName() {
        return SAMPLING_LIKELIHOOD;
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        XMLObject cxo = xo.getChild(MODEL);
        DemographicModel demoModel = (DemographicModel) cxo.getChild(DemographicModel.class);

        cxo = xo.getChild(BETAS);
        Parameter betas = (Parameter) cxo.getChild(Parameter.class);

        cxo = xo.getChild(POPULATION_TREE); // May need to adapt to multiple trees, a la CoalescentLikelihoodParser
        TreeModel tree = (TreeModel) cxo.getChild(TreeModel.class);

        double[] epochWidths = null;
        if (xo.hasChildNamed(EPOCH_WIDTHS)) {
            cxo = xo.getChild(EPOCH_WIDTHS);
            epochWidths = cxo.getDoubleArrayAttribute(WIDTHS);
        }

        MatrixParameter covariates = null;
        if (xo.hasChildNamed(COVARIATES)) {
            cxo = xo.getChild(COVARIATES);
            covariates = (MatrixParameter) cxo.getChild(MatrixParameter.class);
        }

        MatrixParameter powerCovariates = null;
        if (xo.hasChildNamed(POWER_COVARIATES)) {
            cxo = xo.getChild(POWER_COVARIATES);
            powerCovariates = (MatrixParameter) cxo.getChild(MatrixParameter.class);
        }

        Parameter powerBetas = null;
        if (xo.hasChildNamed(POWER_BETAS)) {
            cxo = xo.getChild(POWER_BETAS);
            powerBetas = (Parameter) cxo.getChild(Parameter.class);
        }

        return new BNPRSamplingLikelihood(tree, betas, demoModel, epochWidths, covariates, powerCovariates, powerBetas);
    }

    @Override
    public String getParserDescription() {
        return "This element represents the likelihood of the sampling times given the demographic function.";
    }

    @Override
    public Class getReturnType() {
        return BNPRSamplingLikelihood.class;
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(MODEL, new XMLSyntaxRule[]{
                    new ElementRule(DemographicModel.class)
            }, "The demographic model which describes the effective population size over time"),

            new ElementRule(BETAS, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }, "The log-linear coefficients of effective population size, used to calculate sampling intensity"),

            new ElementRule(POPULATION_TREE, new XMLSyntaxRule[]{
                    new ElementRule(TreeModel.class)
            }, "Tree/sampling times to compute likelihood for"),

            new ElementRule(EPOCH_WIDTHS,
                    new XMLSyntaxRule[]{AttributeRule.newDoubleArrayRule(WIDTHS)}),

            new ElementRule(COVARIATES, new XMLSyntaxRule[]{
                    new ElementRule(MatrixParameter.class)
            }, "Matrix parameter specifying covariate values at latent points.", true),

            new ElementRule(POWER_COVARIATES, new XMLSyntaxRule[]{
                    new ElementRule(MatrixParameter.class)
            }, "Matrix parameter specifying power covariate values at latent points.", true),

            new ElementRule(POWER_BETAS, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }, "The coefficients of the power covariates, used to calculate sampling intensity", true)
    };
}
