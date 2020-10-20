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

import dr.inference.distribution.*;
import dr.inference.model.Parameter;
import dr.inference.operators.repeatedMeasures.GammaGibbsProvider;
import dr.math.MathUtils;
import dr.math.distributions.Distribution;
import dr.math.distributions.GammaDistribution;
import dr.math.matrixAlgebra.Vector;
import dr.xml.*;

/**
 * @author Marc A. Suchard
 * @author Philippe Lemey
 */
public class NormalGammaPrecisionGibbsOperator extends SimpleMCMCOperator implements GibbsOperator, Reportable {

    public static final String OPERATOR_NAME = "normalGammaPrecisionGibbsOperator";
    public static final String LIKELIHOOD = "likelihood";
    public static final String PRIOR = "prior";
    private static final String WORKING = "workingDistribution";

    public NormalGammaPrecisionGibbsOperator(GammaGibbsProvider gammaGibbsProvider, GammaStatisticsProvider prior,
                                             double weight) {
        this(gammaGibbsProvider, prior, null, weight);
    }

    public NormalGammaPrecisionGibbsOperator(GammaGibbsProvider gammaGibbsProvider,
                                             GammaStatisticsProvider prior, GammaStatisticsProvider working,
                                             double weight) {
        this.gammaGibbsProvider = gammaGibbsProvider;
        this.precisionParameter = gammaGibbsProvider.getPrecisionParameter();

        this.prior = prior;
        this.working = working;

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

    @Override
    public String getReport() {
        int dimTrait = precisionParameter.getDimension();
        double[] obsCounts = new double[dimTrait];
        double[] sumSquaredErrors = new double[dimTrait];

        gammaGibbsProvider.drawValues();

        for (int i = 0; i < dimTrait; i++) {
            final GammaGibbsProvider.SufficientStatistics statistics = gammaGibbsProvider.getSufficientStatistics(i);
            obsCounts[i] = statistics.observationCount;
            sumSquaredErrors[i] = statistics.sumOfSquaredErrors;
        }

        StringBuilder sb = new StringBuilder(OPERATOR_NAME + " report:\n");
        sb.append("Observation counts:\t");
        sb.append(new Vector(obsCounts));
        sb.append("\n");
        sb.append("Sum of squared errors:\t");
        sb.append(new Vector(sumSquaredErrors));
        return sb.toString();
    }

    static class GammaParametrization implements GammaStatisticsProvider {
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

        GammaParametrization(Distribution distribution) {
            this(distribution.mean(), distribution.variance());
        }

        double getRate() {
            return rate;
        }

        double getShape() {
            return shape;
        }

        @Override
        public double getShape(int dim) {
            return getShape();
        }

        @Override
        public double getRate(int dim) {
            return getRate();
        }
    }

    private double weigh(double working, double prior) {
        return (1.0 - pathParameter) * working + pathParameter * prior;
    }

    public double doOperation() {

        gammaGibbsProvider.drawValues();

        for (int dim = 0; dim < precisionParameter.getDimension(); ++dim) {

            final GammaGibbsProvider.SufficientStatistics statistics = gammaGibbsProvider.getSufficientStatistics(dim);

            double shape = pathParameter * statistics.observationCount / 2;
            double rate = pathParameter * statistics.sumOfSquaredErrors / 2;

            if (working == null) {

                shape += prior.getShape(dim);
                rate += prior.getRate(dim);

            } else {

                shape += weigh(prior.getShape(dim), prior.getShape(dim)); //TODO: shouldn't these include the working?
                rate += weigh(prior.getRate(dim), prior.getShape(dim));

            }

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

        private void checkGammaDistribution(DistributionLikelihood distribution) throws XMLParseException {
            if (!((distribution.getDistribution() instanceof GammaDistribution) ||
                    (distribution.getDistribution() instanceof GammaDistributionModel))) {
                throw new XMLParseException("Gibbs operator assumes normal-gamma model");
            }
        }

        private GammaStatisticsProvider getGammaStatisticsProvider(Object obj) throws XMLParseException {
            final GammaStatisticsProvider gammaStats;
            if (obj instanceof DistributionLikelihood) {
                DistributionLikelihood priorLike = (DistributionLikelihood) obj;
                checkGammaDistribution(priorLike);
                gammaStats = new GammaParametrization(priorLike.getDistribution());
            } else if (obj instanceof GammaStatisticsProvider) {
                gammaStats = (GammaStatisticsProvider) obj;
            } else {
                throw new XMLParseException("Prior must be gamma");
            }

            return gammaStats;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            final double weight = xo.getDoubleAttribute(WEIGHT);
            final Object prior = xo.getElementFirstChild(PRIOR);
            GammaStatisticsProvider priorDistribution = getGammaStatisticsProvider(prior);


            final Object working = xo.hasChildNamed(WORKING) ? xo.getElementFirstChild(WORKING) : null;

            GammaStatisticsProvider workingDistribution = null;
            if (working != null) {
                workingDistribution = getGammaStatisticsProvider(working);
            }

            final GammaGibbsProvider gammaGibbsProvider;

            if (xo.hasChildNamed(LIKELIHOOD)) {

                DistributionLikelihood likelihood = (DistributionLikelihood) xo.getElementFirstChild(LIKELIHOOD);

                if (!((likelihood.getDistribution() instanceof NormalDistributionModel) ||
                        (likelihood.getDistribution() instanceof LogNormalDistributionModel)
                )) {
                    throw new XMLParseException("Gibbs operator assumes normal-gamma model");
                }

                gammaGibbsProvider = new GammaGibbsProvider.Default(likelihood);

            } else {

                gammaGibbsProvider = (GammaGibbsProvider) xo.getChild(GammaGibbsProvider.class);
            }

            return new NormalGammaPrecisionGibbsOperator(gammaGibbsProvider,
                    priorDistribution, workingDistribution,
                    weight);
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
                                        new XORRule(
                                                new ElementRule(DistributionLikelihood.class),
                                                new ElementRule(GammaStatisticsProvider.class)
                                        )
                                }),

                        new ElementRule(GammaGibbsProvider.class)

                ),
                new ElementRule(PRIOR,
                        new XMLSyntaxRule[]{
                                new ElementRule(DistributionLikelihood.class)
                        }),
                new ElementRule(WORKING,
                        new XMLSyntaxRule[]{
                                new ElementRule(DistributionLikelihood.class)
                        }, true),
        };
    };

    private final GammaGibbsProvider gammaGibbsProvider;
    private final Parameter precisionParameter;

    private final GammaStatisticsProvider prior;
    private final GammaStatisticsProvider working;

    private double pathParameter = 1.0;
}
