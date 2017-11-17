/*
 * PriorParsers.java
 *
 * Copyright (c) 2002-2016 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.inferencexml.distribution;

import dr.inference.distribution.CauchyDistribution;
import dr.inference.distribution.DistributionLikelihood;
import dr.inference.distribution.MultivariateDistributionLikelihood;
import dr.inference.model.Likelihood;
import dr.inference.model.Statistic;
import dr.math.distributions.*;
import dr.xml.*;
import org.apache.commons.math.distribution.CauchyDistributionImpl;

import java.util.ArrayList;
import java.util.List;

/**
 */
public class PriorParsers {
    public final static boolean DEBUG = false;

    public static final String UNIFORM_PRIOR = "uniformPrior";
    public static final String EXPONENTIAL_PRIOR = "exponentialPrior";
    public static final String POISSON_PRIOR = "poissonPrior";
    public static final String NEGATIVE_BINOMIAL_PRIOR = "negativeBinomialPrior";
    public static final String DISCRETE_UNIFORM_PRIOR = "discreteUniformPrior";
    public static final String NORMAL_PRIOR = "normalPrior";
    public static final String LOG_NORMAL_PRIOR = "logNormalPrior";
    public static final String HALF_NORMAL_PRIOR = "halfNormalPrior";
    public static final String GAMMA_PRIOR = "gammaPrior";
    public static final String INVGAMMA_PRIOR = "invgammaPrior";
    public static final String INVGAMMA_PRIOR_CORRECT = "inverseGammaPrior";
    public static final String LAPLACE_PRIOR = "laplacePrior";
    public static final String BETA_PRIOR = "betaPrior";
    public static final String CAUCHY_PRIOR = "cauchyPrior";
    public static final String UPPER = "upper";
    public static final String LOWER = "lower";
    public static final String MEAN = "mean";
    public static final String MEDIAN = "median";
    public static final String STDEV = "stdev";
    public static final String MEAN_IN_REAL_SPACE = "meanInRealSpace";
    public static final String MU = "mu";
    public static final String SIGMA = "sigma";
    public static final String SHAPE = "shape";
    public static final String SHAPEB = "shapeB";
    public static final String SCALE = "scale";
    public static final String DF = "df";
    public static final String OFFSET = "offset";
    public static final String UNINFORMATIVE = "uninformative";
    public static final String HALF_T_PRIOR = "halfTPrior";
    public static final String DIRICHLET_PRIOR = "dirichletPrior";
    public static final String ALPHA = "alpha";
    public static final String COUNTS = "counts";
    public static final String SUMS_TO = "sumsTo";


    /**
     * A special parser that reads a convenient short form of priors on parameters.
     */
    public static XMLObjectParser UNIFORM_PRIOR_PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return UNIFORM_PRIOR;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            double lower = xo.getDoubleAttribute(LOWER);
            double upper = xo.getDoubleAttribute(UPPER);

            if (lower == Double.NEGATIVE_INFINITY || upper == Double.POSITIVE_INFINITY)
                throw new XMLParseException("Uniform prior " + xo.getName() + " cannot take a bound at infinity, " +
                        "because it returns 1/(high-low) = 1/inf");

            DistributionLikelihood likelihood = new DistributionLikelihood(new UniformDistribution(lower, upper));
            if (DEBUG) {
                System.out.println("Uniform prior: " + xo.getChildCount());
            }
            for (int j = 0; j < xo.getChildCount(); j++) {
                if (DEBUG) {
                    System.out.println(xo.getChild(j));
                }
                if (xo.getChild(j) instanceof Statistic) {
                    if (DEBUG) {
                        //System.out.println((Statistic) xo.getChild(j));
                        Statistic test = (Statistic) xo.getChild(j);
                        System.out.println(test.getDimension());
                        for (int i = 0; i < test.getDimension(); i++) {
                            System.out.println("  " + test.getDimensionName(i) + " - " + test.getStatisticValue(i));
                        }
                        System.out.println(test.getClass());
                    }
                    likelihood.addData((Statistic) xo.getChild(j));
                    if (DEBUG) {
                        likelihood.calculateLogLikelihood();
                        System.out.println("likelihood: " + likelihood.getLogLikelihood());
                    }
                } else {
                    throw new XMLParseException("illegal element in " + xo.getName() + " element");
                }
            }
            return likelihood;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                AttributeRule.newDoubleRule(LOWER),
                AttributeRule.newDoubleRule(UPPER),
                new ElementRule(Statistic.class, 1, Integer.MAX_VALUE)
        };

        public String getParserDescription() {
            return "Calculates the prior probability of some data under a given uniform distribution.";
        }

        public Class getReturnType() {
            return Likelihood.class;
        }
    };

    /**
     * A special parser that reads a convenient short form of priors on parameters.
     */
    public static XMLObjectParser EXPONENTIAL_PRIOR_PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return EXPONENTIAL_PRIOR;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            double scale;

            if (xo.hasAttribute(SCALE)) {
                scale = xo.getDoubleAttribute(SCALE);
            } else {
                scale = xo.getDoubleAttribute(MEAN);
            }
            final double offset = xo.hasAttribute(OFFSET) ? xo.getDoubleAttribute(OFFSET) : 0.0;

            DistributionLikelihood likelihood = new DistributionLikelihood(new ExponentialDistribution(1.0 / scale), offset);
            if (DEBUG) {
                System.out.println("Exponential prior: " + xo.getChildCount());
            }
            for (int j = 0; j < xo.getChildCount(); j++) {
                if (DEBUG) {
                    System.out.println(xo.getChild(j));
                }
                if (xo.getChild(j) instanceof Statistic) {
                    if (DEBUG) {
                        System.out.println("scale: " + scale);
                        System.out.println("offset: " + offset);
                        //System.out.println((Statistic) xo.getChild(j));
                        Statistic test = (Statistic) xo.getChild(j);
                        System.out.println(test.getDimension());
                        for (int i = 0; i < test.getDimension(); i++) {
                            System.out.println("  " + test.getDimensionName(i) + " - " + test.getStatisticValue(i));
                        }
                        System.out.println(test.getClass());
                    }
                    likelihood.addData((Statistic) xo.getChild(j));
                    if (DEBUG) {
                        likelihood.makeDirty();
                        likelihood.calculateLogLikelihood();
                        System.out.println("likelihood: " + likelihood.getLogLikelihood());
                    }
                } else {
                    throw new XMLParseException("illegal element in " + xo.getName() + " element");
                }
            }

            return likelihood;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                new XORRule(
                        AttributeRule.newDoubleRule(SCALE),
                        AttributeRule.newDoubleRule(MEAN)
                ),
                AttributeRule.newDoubleRule(OFFSET, true),
                new ElementRule(Statistic.class, 1, Integer.MAX_VALUE)
        };

        public String getParserDescription() {
            return "Calculates the prior probability of some data under a given exponential distribution.";
        }

        public Class getReturnType() {
            return Likelihood.class;
        }
    };

    /**
     * A special parser that reads a convenient short form of priors on parameters.
     */
    public static XMLObjectParser POISSON_PRIOR_PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return POISSON_PRIOR;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            double mean = xo.getDoubleAttribute(MEAN);
            double offset = xo.getDoubleAttribute(OFFSET);

            DistributionLikelihood likelihood = new DistributionLikelihood(new PoissonDistribution(mean), offset);
            for (int j = 0; j < xo.getChildCount(); j++) {
                if (xo.getChild(j) instanceof Statistic) {
                    likelihood.addData((Statistic) xo.getChild(j));
                } else {
                    throw new XMLParseException("illegal element in " + xo.getName() + " element");
                }
            }

            return likelihood;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                AttributeRule.newDoubleRule(MEAN),
                AttributeRule.newDoubleRule(OFFSET),
                new ElementRule(Statistic.class, 1, Integer.MAX_VALUE)
        };

        public String getParserDescription() {
            return "Calculates the prior probability of some data under a given poisson distribution.";
        }

        public Class getReturnType() {
            return Likelihood.class;
        }
    };

    /**
     * A special parser that reads a convenient short form of priors on parameters.
     */
    public static XMLObjectParser NEGATIVE_BINOMIAL_PRIOR_PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return NEGATIVE_BINOMIAL_PRIOR;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            double mean = xo.getDoubleAttribute(MEAN);
            double stdev = xo.getDoubleAttribute(STDEV);

            DistributionLikelihood likelihood = new DistributionLikelihood(new NegativeBinomialDistribution(mean, stdev));
            for (int j = 0; j < xo.getChildCount(); j++) {
                if (xo.getChild(j) instanceof Statistic) {
                    likelihood.addData((Statistic) xo.getChild(j));
                } else {
                    throw new XMLParseException("illegal element in " + xo.getName() + " element");
                }
            }

            return likelihood;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                AttributeRule.newDoubleRule(MEAN),
                AttributeRule.newDoubleRule(STDEV),
                new ElementRule(Statistic.class, 1, Integer.MAX_VALUE)
        };

        public String getParserDescription() {
            return "Calculates the prior probability of some data under a given negative binomial distribution.";
        }

        public Class getReturnType() {
            return Likelihood.class;
        }
    };

    /**
     * A special parser that reads a convenient short form of priors on parameters.
     */
    public static XMLObjectParser DISCRETE_UNIFORM_PRIOR_PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return DISCRETE_UNIFORM_PRIOR;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            double upper = xo.getDoubleAttribute(UPPER);
            double lower = xo.getDoubleAttribute(LOWER);

            DistributionLikelihood likelihood = new DistributionLikelihood(new DiscreteUniformDistribution(lower, upper));
            for (int j = 0; j < xo.getChildCount(); j++) {
                if (xo.getChild(j) instanceof Statistic) {
                    likelihood.addData((Statistic) xo.getChild(j));
                } else {
                    throw new XMLParseException("illegal element in " + xo.getName() + " element");
                }
            }

            return likelihood;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                AttributeRule.newDoubleRule(LOWER),
                AttributeRule.newDoubleRule(UPPER),
                new ElementRule(Statistic.class, 1, Integer.MAX_VALUE)
        };

        public String getParserDescription() {
            return "Calculates the prior probability of some data under a given discrete uniform distribution.";
        }

        public Class getReturnType() {
            return Likelihood.class;
        }
    };



    /**
     * A special parser that reads a convenient short form of priors on parameters.
     */
    public static XMLObjectParser HALF_T_PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return HALF_T_PRIOR;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            double scale = xo.getDoubleAttribute(SCALE);
            double df = xo.getDoubleAttribute(DF);

            DistributionLikelihood likelihood = new DistributionLikelihood(new HalfTDistribution(scale, df));
            for (int j = 0; j < xo.getChildCount(); j++) {
                if (xo.getChild(j) instanceof Statistic) {
                    likelihood.addData((Statistic) xo.getChild(j));
                } else {
                    throw new XMLParseException("illegal element in " + xo.getName() + " element");
                }
            }

            return likelihood;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                AttributeRule.newDoubleRule(SCALE),
                AttributeRule.newDoubleRule(DF),
                new ElementRule(Statistic.class, 1, Integer.MAX_VALUE)
        };

        public String getParserDescription() {
            return "Calculates the prior probability of some data under a given half-T distribution.";
        }

        public Class getReturnType() {
            return Likelihood.class;
        }
    };

    /**
     * A special parser that reads a convenient short form of priors on parameters.
     */
    public static XMLObjectParser NORMAL_PRIOR_PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return NORMAL_PRIOR;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            double mean = xo.getDoubleAttribute(MEAN);
            double stdev = xo.getDoubleAttribute(STDEV);

            DistributionLikelihood likelihood = new DistributionLikelihood(new NormalDistribution(mean, stdev));
            for (int j = 0; j < xo.getChildCount(); j++) {
                if (xo.getChild(j) instanceof Statistic) {
                    likelihood.addData((Statistic) xo.getChild(j));
                } else {
                    throw new XMLParseException("illegal element in " + xo.getName() + " element");
                }
            }

            return likelihood;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                AttributeRule.newDoubleRule(MEAN),
                AttributeRule.newDoubleRule(STDEV),
                new ElementRule(Statistic.class, 1, Integer.MAX_VALUE)
        };

        public String getParserDescription() {
            return "Calculates the prior probability of some data under a given normal distribution.";
        }

        public Class getReturnType() {
            return Likelihood.class;
        }
    };

    public static XMLObjectParser HALF_NORMAL_PRIOR_PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return HALF_NORMAL_PRIOR;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            double mean = xo.getDoubleAttribute(MEAN);
            double stdev = xo.getDoubleAttribute(STDEV);

            DistributionLikelihood likelihood = new DistributionLikelihood(new HalfNormalDistribution(mean, stdev));
            for (int j = 0; j < xo.getChildCount(); j++) {
                if (xo.getChild(j) instanceof Statistic) {
                    likelihood.addData((Statistic) xo.getChild(j));
                } else {
                    throw new XMLParseException("illegal element in " + xo.getName() + " element");
                }
            }

            return likelihood;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                AttributeRule.newDoubleRule(MEAN),
                AttributeRule.newDoubleRule(STDEV),
                new ElementRule(Statistic.class, 1, Integer.MAX_VALUE)
        };

        public String getParserDescription() {
            return "Calculates the prior probability of some data under a given half-normal distribution.";
        }

        public Class getReturnType() {
            return Likelihood.class;
        }
    };


    /**
     * A special parser that reads a convenient short form of priors on parameters.
     * <p/>
     * If X ~ logNormal, then log(X) ~ Normal.
     * <br>
     * <br>
     * If meanInRealSpace=false, <code>mean</code> specifies the mean of log(X) and
     * <code>stdev</code> specifies the standard deviation of log(X).
     * <br>
     * <br>
     * If meanInRealSpace=true, <code>mean</code> specifies the mean of X, but <code>
     * stdev</code> specifies the standard deviation of log(X).
     * <br>
     */
    public static XMLObjectParser LOG_NORMAL_PRIOR_PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return LOG_NORMAL_PRIOR;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            double mu;
            double sigma;

            final boolean meanInRealSpace = xo.getAttribute(MEAN_IN_REAL_SPACE, false);

            if (xo.hasAttribute(MEAN_IN_REAL_SPACE)) {
                // if the meanInRealSpace attribute is given then respect the old
                // parser choices (i.e., specify mean and stdev but change their
                // meaning to mu and sigma).

                final double mean = xo.getDoubleAttribute(MEAN);
                final double stdev = xo.getDoubleAttribute(STDEV);

                if (meanInRealSpace) {
                    mu = Math.log(mean / Math.sqrt(1 + (stdev * stdev) / (mean * mean)));
                } else {
                    mu = mean;
                }

                // in the previous version, the stdev (or 'Log(stdev)') parameter is sigma.
                sigma = stdev;
            } else {
                // decide parameterization by whether mean, stdev or mu, sigma are given.
                if (xo.hasAttribute(MEAN)) {
                    final double mean = xo.getDoubleAttribute(MEAN);
                    if (!xo.hasAttribute(STDEV)) {
                        throw new XMLParseException("Lognormal should be specified either with mean and stdev or mu and sigma.");
                    }
                    final double stdev = xo.getDoubleAttribute(STDEV);
                    if (mean <= 0) {
                        throw new XMLParseException("If specified with a mean in real space, the value should be positive.");
                    }
                    mu = Math.log(mean / Math.sqrt(1 + (stdev * stdev) / (mean * mean)));
                    sigma = Math.sqrt(Math.log(1 + (stdev * stdev) / (mean * mean)));
                } else {
                    if (meanInRealSpace) {
                        throw new XMLParseException("Lognormal with 'meanInRealSpace' should be specified with mean and stdev.");
                    }
                    mu = xo.getDoubleAttribute(MU);
                    if (!xo.hasAttribute(SIGMA)) {
                        throw new XMLParseException("Lognormal should be specified either with mean and stdev or mu and sigma.");
                    }
                    sigma = xo.getDoubleAttribute(SIGMA);
                }
            }
            final double offset = xo.getAttribute(OFFSET, 0.0);

            final DistributionLikelihood likelihood = new DistributionLikelihood(new LogNormalDistribution(mu, sigma), offset);

            for (int j = 0; j < xo.getChildCount(); j++) {
                if (xo.getChild(j) instanceof Statistic) {
                    likelihood.addData((Statistic) xo.getChild(j));
                } else {
                    throw new XMLParseException("illegal element in " + xo.getName() + " element");
                }
            }

            return likelihood;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                new XORRule(
                        AttributeRule.newDoubleRule(MEAN),
                        AttributeRule.newDoubleRule(MU)
                ),
                new XORRule(
                        AttributeRule.newDoubleRule(STDEV),
                        AttributeRule.newDoubleRule(SIGMA)
                ),
                AttributeRule.newDoubleRule(OFFSET, true),
                AttributeRule.newBooleanRule(MEAN_IN_REAL_SPACE, true),
                new ElementRule(Statistic.class, 1, Integer.MAX_VALUE)
        };

        public String getParserDescription() {
            return "Calculates the prior probability of some data under a given lognormal distribution.";
        }

        public Class getReturnType() {
            return Likelihood.class;
        }
    };


    /**
     * A special parser that reads a convenient short form of priors on parameters.
     */
    public static XMLObjectParser GAMMA_PRIOR_PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return GAMMA_PRIOR;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            final double shape = xo.getDoubleAttribute(SHAPE);
            final double scale = xo.getDoubleAttribute(SCALE);
            final double offset = xo.getAttribute(OFFSET, 0.0);

            DistributionLikelihood likelihood = new DistributionLikelihood(new GammaDistribution(shape, scale), offset);
            for (int j = 0; j < xo.getChildCount(); j++) {
                if (xo.getChild(j) instanceof Statistic) {
                    likelihood.addData((Statistic) xo.getChild(j));
                } else {
                    throw new XMLParseException("illegal element in " + xo.getName() + " element");
                }
            }

            return likelihood;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                AttributeRule.newDoubleRule(SHAPE),
                AttributeRule.newDoubleRule(SCALE),
                AttributeRule.newDoubleRule(OFFSET, true),
                // AttributeRule.newBooleanRule(UNINFORMATIVE, true),
                new ElementRule(Statistic.class, 1, Integer.MAX_VALUE)
        };

        public String getParserDescription() {
            return "Calculates the prior probability of some data under a given gamma distribution.";
        }

        public Class getReturnType() {
            return Likelihood.class;
        }
    };

    /**
     * A special parser that reads a convenient short form of priors on parameters.
     */
    public static XMLObjectParser INVGAMMA_PRIOR_PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return INVGAMMA_PRIOR;
        }

        public String[] getParserNames() {
            return new String[]{INVGAMMA_PRIOR, INVGAMMA_PRIOR_CORRECT};
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            final double shape = xo.getDoubleAttribute(SHAPE);
            final double scale = xo.getDoubleAttribute(SCALE);
            final double offset = xo.getAttribute(OFFSET, 0.0);

            DistributionLikelihood likelihood = new DistributionLikelihood(new InverseGammaDistribution(shape, scale), offset);

            for (int j = 0; j < xo.getChildCount(); j++) {
                if (xo.getChild(j) instanceof Statistic) {
                    likelihood.addData((Statistic) xo.getChild(j));
                } else {
                    throw new XMLParseException("illegal element in " + xo.getName() + " element");
                }
            }

            return likelihood;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                AttributeRule.newDoubleRule(SHAPE),
                AttributeRule.newDoubleRule(SCALE),
                AttributeRule.newDoubleRule(OFFSET, true),
                new ElementRule(Statistic.class, 1, Integer.MAX_VALUE)
        };

        public String getParserDescription() {
            return "Calculates the prior probability of some data under a given inverse gamma distribution.";
        }

        public Class getReturnType() {
            return Likelihood.class;
        }
    };

    public static XMLObjectParser LAPLACE_PRIOR_PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return LAPLACE_PRIOR;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            double mean = xo.getDoubleAttribute(MEAN);
            double scale = xo.getDoubleAttribute(SCALE);

            DistributionLikelihood likelihood = new DistributionLikelihood(new LaplaceDistribution(mean, scale));
            for (int j = 0; j < xo.getChildCount(); j++) {
                if (xo.getChild(j) instanceof Statistic) {
                    likelihood.addData((Statistic) xo.getChild(j));
                } else {
                    throw new XMLParseException("illegal element in " + xo.getName() + " element");
                }
            }

            return likelihood;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                AttributeRule.newDoubleRule(MEAN),
                AttributeRule.newDoubleRule(SCALE),
                new ElementRule(Statistic.class, 1, Integer.MAX_VALUE)
        };

        public String getParserDescription() {
            return "Calculates the prior probability of some data under a given laplace distribution.";
        }

        public Class getReturnType() {
            return Likelihood.class;
        }
    };

    /**
     * A special parser that reads a convenient short form of priors on parameters.
     */
    public static XMLObjectParser BETA_PRIOR_PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return BETA_PRIOR;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            final double shape = xo.getDoubleAttribute(SHAPE);
            final double shapeB = xo.getDoubleAttribute(SHAPEB);
            final double offset = xo.getAttribute(OFFSET, 0.0);
            final double scale = xo.getAttribute(SCALE, 1.0);

            DistributionLikelihood likelihood = new DistributionLikelihood(new BetaDistribution(shape, shapeB), offset, scale);
            for (int j = 0; j < xo.getChildCount(); j++) {
                if (xo.getChild(j) instanceof Statistic) {
                    likelihood.addData((Statistic) xo.getChild(j));
                } else {
                    throw new XMLParseException("illegal element in " + xo.getName() + " element");
                }
            }

            return likelihood;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                AttributeRule.newDoubleRule(SHAPE),
                AttributeRule.newDoubleRule(SHAPEB),
                AttributeRule.newDoubleRule(OFFSET, true),
                new ElementRule(Statistic.class, 1, Integer.MAX_VALUE)
        };

        public String getParserDescription() {
            return "Calculates the prior probability of some data under a given beta distribution.";
        }

        public Class getReturnType() {
            return Likelihood.class;
        }
    };

    /**
     * A special parser that reads a convenient short form of priors on parameters.
     */
    public static XMLObjectParser DIRICHLET_PRIOR_PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return DIRICHLET_PRIOR;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            List<Statistic> data = new ArrayList<Statistic>();
            int dim = 0;
            double sum = 0;
            for (int j = 0; j < xo.getChildCount(); j++) {
                if (xo.getChild(j) instanceof Statistic) {
                    Statistic statistic = (Statistic) xo.getChild(j);
                    dim += statistic.getDimension();
                    sum += statistic.getValueSum();
                    data.add(statistic);
                } else {
                    throw new XMLParseException("Illegal element in " + xo.getName() + " element.");
                }
            }

            double[] counts;

            if (xo.hasAttribute(ALPHA)) {
                // alpha is for when the Dirichlet distribution is symmetrical
                double alpha = xo.getDoubleAttribute(ALPHA);
                counts  = new double[dim];
                for (int i = 0; i < dim; i++) {
                    counts[i] = alpha;
                }
            } else {
                counts = xo.getDoubleArrayAttribute(COUNTS);
                if (counts.length != dim) {
                    throw new XMLParseException("Counts attribute in " + xo.getName() + " should have the same dimension as the data.");
                }
            }

            double sumsTo = 1.0;
            if (xo.hasAttribute(SUMS_TO)) {
                sumsTo = xo.getDoubleAttribute(SUMS_TO);
            }
            if (Math.abs(sumsTo - sum) > DirichletDistribution.ACCURACY_THRESHOLD) {
                throw new XMLParseException("The starting values of the data should sum to " + sumsTo + " (actually sums to " + sum + ").");
            }

            MultivariateDistributionLikelihood likelihood = new MultivariateDistributionLikelihood(new DirichletDistribution(counts, sumsTo));

            for (Statistic statistic : data) {
                likelihood.addData(statistic);
            }

            return likelihood;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                new XORRule(
                        AttributeRule.newDoubleArrayRule(ALPHA),
                        AttributeRule.newDoubleArrayRule(COUNTS)
                ),
                AttributeRule.newDoubleRule(SUMS_TO, true),
                new ElementRule(Statistic.class, 1, Integer.MAX_VALUE)
        };

        public String getParserDescription() {
            return "Calculates the prior probability of some data under a Dirichlet distribution.";
        }

        public Class getReturnType() {
            return MultivariateDistributionLikelihood.class;
        }
    };



    public static XMLObjectParser CAUCH_PRIOR_PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return CAUCHY_PRIOR;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            final double median = xo.getDoubleAttribute(MEDIAN);
            final double scale = xo.getDoubleAttribute(SCALE);

            DistributionLikelihood likelihood = new DistributionLikelihood(new CauchyDistribution(median, scale));
            for (int j = 0; j < xo.getChildCount(); j++) {
                if (xo.getChild(j) instanceof Statistic) {
                    likelihood.addData((Statistic) xo.getChild(j));
                } else {
                    throw new XMLParseException("illegal element in " + xo.getName() + " element");
                }
            }

            return likelihood;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                AttributeRule.newDoubleRule(MEDIAN),
                AttributeRule.newDoubleRule(SCALE),
                new ElementRule(Statistic.class, 1, Integer.MAX_VALUE)
        };
        public String getParserDescription() {
            return "Calculates the prior probability of some data under a Cauchy distribution.";
        }

        public Class getReturnType() {
            return DistributionLikelihood.class;
        }
    };

}
