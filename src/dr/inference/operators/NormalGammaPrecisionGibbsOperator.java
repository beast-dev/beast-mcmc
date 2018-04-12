/*
 * NormalGammaPrecisionGibbsOperator.java
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

package dr.inference.operators;

import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.continuous.RepeatedMeasuresTraitDataModel;
import dr.inference.distribution.DistributionLikelihood;
import dr.inference.distribution.GammaDistributionModel;
import dr.inference.distribution.LogNormalDistributionModel;
import dr.inference.distribution.NormalDistributionModel;
import dr.inference.model.Parameter;
import dr.inference.operators.repeatedMeasures.dr.inference.operators.repeatedMeasures.GammaGibbsProvider;
import dr.math.MathUtils;
import dr.math.distributions.Distribution;
import dr.math.distributions.GammaDistribution;
import dr.xml.*;

/**
 * @author Marc A. Suchard
 * @author Philippe Lemey
 */
public class NormalGammaPrecisionGibbsOperator extends SimpleMCMCOperator implements GibbsOperator {

    public static final String OPERATOR_NAME = "normalGammaPrecisionGibbsOperator";
    public static final String LIKELIHOOD = "likelihood";
    private static final String REPEATED_MEASURES = "repeatedMeasures";
    public static final String PRIOR = "prior";
    
    public NormalGammaPrecisionGibbsOperator(GammaGibbsProvider gammaGibbsProvider, Distribution prior,
                                             double weight) {
        this.gammaGibbsProvider = gammaGibbsProvider;
        this.precisionParameter = gammaGibbsProvider.getPrecisionParameter();
        this.prior = prior;

        setWeight(weight);
    }

    /**
     * @return a short descriptive message of the performance of this operator.
     */
    public String getPerformanceSuggestion() {
        return null;
    }

    public String getOperatorName() {
        return OPERATOR_NAME;
    }

    static class GammaParametrization {
        private final double rate;
        private final double shape;

        GammaParametrization(double mean, double variance) {
            if (mean == 0) {
                rate = 0;
                shape = -0.5; // Uninformative prior
            } else {
                rate = mean / variance;
                shape = mean * rate;
            }
        }

        double getRate() { return rate; }
        double getShape() { return shape; }
    }

    public double doOperation() {

        GammaParametrization priorParametrization = new GammaParametrization(
                prior.mean(),
                prior.variance());

        gammaGibbsProvider.drawValues();

        for (int dim = 0; dim < precisionParameter.getDimension(); ++dim) {

            final GammaGibbsProvider.SufficientStatistics statistics = gammaGibbsProvider.getSufficientStatistics(dim);

            final double shape = priorParametrization.getShape() + pathParameter * statistics.observationCount / 2;
            final double rate = priorParametrization.getRate() + pathParameter * statistics.sumOfSquaredErrors / 2;

            final double draw = MathUtils.nextGamma(shape, rate); // Gamma( \alpha + n/2 , \beta + (1/2)*SSE )

            precisionParameter.setParameterValue(dim, draw);
        }

        return 0;
    }

    @Override
    public void setPathParameter(double beta) {
        if (beta < 0.0 || beta > 1.0) {
            throw new IllegalArgumentException("Invalid pathParameter value");
        }

        this.pathParameter = beta;
    }

    /**
     * @return the number of steps the operator performs in one go.
     */
    public int getStepCount() {
        return 1;
    }

    public static dr.xml.XMLObjectParser PARSER = new dr.xml.AbstractXMLObjectParser() {

        public String getParserName() {
            return OPERATOR_NAME;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            double weight = xo.getDoubleAttribute(WEIGHT);

//            DistributionLikelihood prior = (DistributionLikelihood) xo.getElementFirstChild(PRIOR);
//
//            if (!((prior.getDistribution() instanceof GammaDistribution) ||
//                    (prior.getDistribution() instanceof GammaDistributionModel)
//            ) ||
//                    !((likelihood.getDistribution() instanceof NormalDistributionModel) ||
//                            (likelihood.getDistribution() instanceof LogNormalDistributionModel)
//                    ))
//                throw new XMLParseException("Gibbs operator assumes normal-gamma model");

            DistributionLikelihood prior = (DistributionLikelihood) xo.getElementFirstChild(PRIOR);

            if (!((prior.getDistribution() instanceof GammaDistribution) ||
                    (prior.getDistribution() instanceof GammaDistributionModel))) {
                throw new XMLParseException("Gibbs operator assumes normal-gamma model");
            }

            GammaGibbsProvider gammaGibbsProvider;

            if (xo.hasChildNamed(LIKELIHOOD)) {

                DistributionLikelihood likelihood = (DistributionLikelihood) xo.getElementFirstChild(LIKELIHOOD);

                if (!((likelihood.getDistribution() instanceof NormalDistributionModel) ||
                                            (likelihood.getDistribution() instanceof LogNormalDistributionModel)
                                    )) {
                    throw new XMLParseException("Gibbs operator assumes normal-gamma model");
                }

                gammaGibbsProvider = new GammaGibbsProvider.Default(likelihood);

            } else {

                XMLObject cxo = xo.getChild(REPEATED_MEASURES);

                RepeatedMeasuresTraitDataModel dataModel = (RepeatedMeasuresTraitDataModel)
                        cxo.getChild(RepeatedMeasuresTraitDataModel.class);

                TreeDataLikelihood likelihood = (TreeDataLikelihood) cxo.getChild(TreeDataLikelihood.class);

                gammaGibbsProvider = new GammaGibbsProvider.RepeatedMeasuresGibbsProvider(
                        dataModel, likelihood, dataModel.getTraitName());
            }

            return new NormalGammaPrecisionGibbsOperator(gammaGibbsProvider, prior.getDistribution(), weight);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "This element returns a operator on the precision parameter of a normal model with gamma prior.";
        }

        public Class getReturnType() {
            return MCMCOperator.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                AttributeRule.newDoubleRule(WEIGHT),
                new XORRule(
                        new ElementRule(LIKELIHOOD,
                                new XMLSyntaxRule[]{
                                        new ElementRule(DistributionLikelihood.class)
                                }),
                        new ElementRule(REPEATED_MEASURES,
                                new XMLSyntaxRule[]{
                                        new ElementRule(RepeatedMeasuresTraitDataModel.class),
                                        new ElementRule(TreeDataLikelihood.class),
                                })
                ),
                new ElementRule(PRIOR,
                        new XMLSyntaxRule[]{
                                new ElementRule(DistributionLikelihood.class)
                        }),
        };
    };

    private final GammaGibbsProvider gammaGibbsProvider;
    private final Distribution prior;
    private final Parameter precisionParameter;

    private double pathParameter = 1.0;
}
