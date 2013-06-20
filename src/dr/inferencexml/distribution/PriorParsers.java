/*
 * PriorParsers.java
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

package dr.inferencexml.distribution;

import dr.evolution.io.Importer.ImportException;
import dr.evolution.io.TreeTrace;
import dr.evomodel.tree.ConditionalCladeFrequency;
import dr.evomodel.tree.TreeModel;
import dr.inference.distribution.ConditionalCladeProbability;
import dr.inference.distribution.DistributionLikelihood;
import dr.inference.distribution.MultivariateDistributionLikelihood;
import dr.inference.model.ConstantLikelihood;
import dr.inference.model.Likelihood;
import dr.inference.model.Statistic;
import dr.inference.trace.LogFileTraces;
import dr.inference.trace.TraceException;
import dr.math.distributions.*;
import dr.util.FileHelpers;
import dr.xml.*;
import jebl.math.Binomial;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.Reader;

/**
 */
public class PriorParsers {
    public static final String UNIFORM_PRIOR = "uniformPrior";
    public static final String EXPONENTIAL_PRIOR = "exponentialPrior";
    public static final String POISSON_PRIOR = "poissonPrior";
    public static final String NORMAL_PRIOR = "normalPrior";
    public static final String NORMAL_REFERENCE_PRIOR = "normalReferencePrior";
    public static final String CONDITIONAL_CLADE_REFERENCE_PRIOR = "conditionalCladeProbability";
    public static final String COALESCENT_HEIGHTS_REFERENCE_PRIOR = "coalescentHeightsReferencePrior";
    public final static String BURNIN = "burnin";
    public final static String EPSILON = "epsilon";
    public static final String LOG_TRANSFORMED_NORMAL_REFERENCE_PRIOR = "logTransformedNormalReferencePrior";
    public static final String LOG_NORMAL_PRIOR = "logNormalPrior";
    public static final String GAMMA_PRIOR = "gammaPrior";
    public static final String GAMMA_REFERENCE_PRIOR = "gammaReferencePrior";
    public static final String INVGAMMA_PRIOR = "invgammaPrior";
    public static final String INVGAMMA_PRIOR_CORRECT = "inverseGammaPrior";
    public static final String LAPLACE_PRIOR = "laplacePrior";
    public static final String BETA_PRIOR = "betaPrior";
    public static final String PARAMETER_COLUMN = "parameterColumn";
    public static final String UPPER = "upper";
    public static final String LOWER = "lower";
    public static final String MEAN = "mean";
    public static final String MEAN_IN_REAL_SPACE = "meanInRealSpace";
    public static final String STDEV = "stdev";
    public static final String SHAPE = "shape";
    public static final String SHAPEB = "shapeB";
    public static final String SCALE = "scale";
    public static final String DF = "df";
    public static final String OFFSET = "offset";
    public static final String UNINFORMATIVE = "uninformative";
    public static final String HALF_T_PRIOR = "halfTPrior";

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
            System.out.println("Uniform prior: " + xo.getChildCount());
            for (int j = 0; j < xo.getChildCount(); j++) {
                System.out.println(xo.getChild(j));
                if (xo.getChild(j) instanceof Statistic) {
                    //System.out.println((Statistic) xo.getChild(j));
                    Statistic test = (Statistic) xo.getChild(j);
                    System.out.println(test.getDimension());
                    for (int i = 0; i < test.getDimension(); i++) {
                        System.out.println("  " + test.getDimensionName(i) + " - " + test.getStatisticValue(i));
                    }
                    System.out.println(test.getClass());
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
            System.out.println("Exponential prior: " + xo.getChildCount());
            for (int j = 0; j < xo.getChildCount(); j++) {
                System.out.println(xo.getChild(j));
                if (xo.getChild(j) instanceof Statistic) {
                    //System.out.println((Statistic) xo.getChild(j));
                    Statistic test = (Statistic) xo.getChild(j);
                    System.out.println(test.getDimension());
                    for (int i = 0; i < test.getDimension(); i++) {
                        System.out.println("  " + test.getDimensionName(i) + " - " + test.getStatisticValue(i));
                    }
                    System.out.println(test.getClass());
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

    /**
     * A special parser that reads a convenient short form of reference priors on parameters.
     */
    public static XMLObjectParser GAMMA_REFERENCE_PRIOR_PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return GAMMA_REFERENCE_PRIOR;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            String fileName = xo.getStringAttribute(FileHelpers.FILE_NAME);

            try {

                File file = new File(fileName);
                String parent = file.getParent();

                if (!file.isAbsolute()) {
                    parent = System.getProperty("user.dir");
                }
                file = new File(parent, fileName);
                fileName = file.getAbsolutePath();

                String parameterName = xo.getStringAttribute(PARAMETER_COLUMN);

                LogFileTraces traces = new LogFileTraces(fileName, file);
                traces.loadTraces();
                int maxState = traces.getMaxState();

                // leaving the burnin attribute off will result in 10% being used
                int burnin = xo.getAttribute("burnin", maxState / 10);
                if (burnin < 0 || burnin >= maxState) {
                    burnin = maxState / 10;
                    System.out.println("WARNING: Burn-in larger than total number of states - using 10%");
                }
                traces.setBurnIn(burnin);

                int traceIndexParameter = -1;
                for (int i = 0; i < traces.getTraceCount(); i++) {
                    String traceName = traces.getTraceName(i);
                    if (traceName.trim().equals(parameterName)) {
                        traceIndexParameter = i;
                    }
                }

                if (traceIndexParameter == -1) {
                    throw new XMLParseException("Column '" + parameterName + "' can not be found for " + getParserName() + " element.");
                }

                Double[] parameterSamples = new Double[traces.getStateCount()];

                DistributionLikelihood likelihood = new DistributionLikelihood(new GammaKDEDistribution((Double[]) traces.getValues(traceIndexParameter).toArray(parameterSamples)));
                for (int j = 0; j < xo.getChildCount(); j++) {
                    if (xo.getChild(j) instanceof Statistic) {
                        likelihood.addData((Statistic) xo.getChild(j));
                    } else {
                        throw new XMLParseException("illegal element in " + xo.getName() + " element");
                    }
                }

                return likelihood;

            } catch (FileNotFoundException fnfe) {
                throw new XMLParseException("File '" + fileName + "' can not be opened for " + getParserName() + " element.");
            } catch (java.io.IOException ioe) {
                throw new XMLParseException(ioe.getMessage());
            } catch (TraceException e) {
                throw new XMLParseException(e.getMessage());
            }
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                AttributeRule.newStringRule("fileName"),
                AttributeRule.newStringRule("parameterColumn"),
                AttributeRule.newIntegerRule("burnin"),
                new ElementRule(Statistic.class, 1, Integer.MAX_VALUE)
        };

        public String getParserDescription() {
            return "Calculates the reference prior probability of some data under a given normal distribution.";
        }

        public Class getReturnType() {
            return Likelihood.class;
        }
    };

    /**
     * A special parser that reads a convenient short form of reference priors on parameters.
     */
    public static XMLObjectParser LOG_TRANSFORMED_NORMAL_REFERENCE_PRIOR_PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return LOG_TRANSFORMED_NORMAL_REFERENCE_PRIOR;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            String fileName = xo.getStringAttribute(FileHelpers.FILE_NAME);

            try {

                File file = new File(fileName);
                String parent = file.getParent();

                if (!file.isAbsolute()) {
                    parent = System.getProperty("user.dir");
                }
                file = new File(parent, fileName);
                fileName = file.getAbsolutePath();

                String parameterName = xo.getStringAttribute(PARAMETER_COLUMN);

                LogFileTraces traces = new LogFileTraces(fileName, file);
                traces.loadTraces();
                int maxState = traces.getMaxState();

                // leaving the burnin attribute off will result in 10% being used
                int burnin = xo.getAttribute("burnin", maxState / 10);
                if (burnin < 0 || burnin >= maxState) {
                    burnin = maxState / 10;
                    System.out.println("WARNING: Burn-in larger than total number of states - using 10%");
                }
                traces.setBurnIn(burnin);

                int traceIndexParameter = -1;
                for (int i = 0; i < traces.getTraceCount(); i++) {
                    String traceName = traces.getTraceName(i);
                    if (traceName.trim().equals(parameterName)) {
                        traceIndexParameter = i;
                    }
                }

                if (traceIndexParameter == -1) {
                    throw new XMLParseException("Column '" + parameterName + "' can not be found for " + getParserName() + " element.");
                }

                Double[] parameterSamples = new Double[traces.getStateCount()];
                traces.getValues(traceIndexParameter).toArray(parameterSamples);

                DistributionLikelihood likelihood = new DistributionLikelihood(new LogTransformedNormalKDEDistribution(parameterSamples));
                for (int j = 0; j < xo.getChildCount(); j++) {
                    if (xo.getChild(j) instanceof Statistic) {
                        likelihood.addData((Statistic) xo.getChild(j));
                    } else {
                        throw new XMLParseException("illegal element in " + xo.getName() + " element");
                    }
                }

                return likelihood;

            } catch (FileNotFoundException fnfe) {
                throw new XMLParseException("File '" + fileName + "' can not be opened for " + getParserName() + " element.");
            } catch (java.io.IOException ioe) {
                throw new XMLParseException(ioe.getMessage());
            } catch (TraceException e) {
                throw new XMLParseException(e.getMessage());
            }

        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                AttributeRule.newStringRule("fileName"),
                AttributeRule.newStringRule("parameterColumn"),
                AttributeRule.newIntegerRule("burnin"),
                new ElementRule(Statistic.class, 1, Integer.MAX_VALUE)
        };

        public String getParserDescription() {
            return "Calculates the reference prior probability of some data under log transformed normal distribution.";
        }

        public Class getReturnType() {
            return Likelihood.class;
        }

    };

    /**
     * A special parser that reads a convenient short form of reference priors on parameters.
     */
    public static XMLObjectParser NORMAL_REFERENCE_PRIOR_PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return NORMAL_REFERENCE_PRIOR;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            String fileName = xo.getStringAttribute(FileHelpers.FILE_NAME);

            try {

                File file = new File(fileName);
                String parent = file.getParent();

                if (!file.isAbsolute()) {
                    parent = System.getProperty("user.dir");
                }
                file = new File(parent, fileName);
                fileName = file.getAbsolutePath();

                String parameterName = xo.getStringAttribute(PARAMETER_COLUMN);

                LogFileTraces traces = new LogFileTraces(fileName, file);
                traces.loadTraces();
                int maxState = traces.getMaxState();

                // leaving the burnin attribute off will result in 10% being used
                int burnin = xo.getAttribute("burnin", maxState / 10);
                if (burnin < 0 || burnin >= maxState) {
                    burnin = maxState / 10;
                    System.out.println("WARNING: Burn-in larger than total number of states - using 10%");
                }
                traces.setBurnIn(burnin);

                int traceIndexParameter = -1;
                for (int i = 0; i < traces.getTraceCount(); i++) {
                    String traceName = traces.getTraceName(i);
                    if (traceName.trim().equals(parameterName)) {
                        traceIndexParameter = i;
                    }
                }

                if (traceIndexParameter == -1) {
                    throw new XMLParseException("Column '" + parameterName + "' can not be found for " + getParserName() + " element.");
                }

                Double[] parameterSamples = new Double[traces.getStateCount()];

                DistributionLikelihood likelihood = new DistributionLikelihood(new NormalKDEDistribution((Double[]) traces.getValues(traceIndexParameter).toArray(parameterSamples)));
                for (int j = 0; j < xo.getChildCount(); j++) {
                    if (xo.getChild(j) instanceof Statistic) {
                        likelihood.addData((Statistic) xo.getChild(j));
                    } else {
                        throw new XMLParseException("illegal element in " + xo.getName() + " element");
                    }
                }

                return likelihood;

            } catch (FileNotFoundException fnfe) {
                throw new XMLParseException("File '" + fileName + "' can not be opened for " + getParserName() + " element.");
            } catch (java.io.IOException ioe) {
                throw new XMLParseException(ioe.getMessage());
            } catch (TraceException e) {
                throw new XMLParseException(e.getMessage());
            }
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                AttributeRule.newStringRule("fileName"),
                AttributeRule.newStringRule("parameterColumn"),
                AttributeRule.newIntegerRule("burnin"),
                new ElementRule(Statistic.class, 1, Integer.MAX_VALUE)
        };

        public String getParserDescription() {
            return "Calculates the reference prior probability of some data under a given normal distribution.";
        }

        public Class getReturnType() {
            return Likelihood.class;
        }
    };

    /**
     * A special parser that reads coalescent heights and builds reference priors for them.
     */
    public static XMLObjectParser COALESCENT_HEIGHTS_REFERENCE_PRIOR_PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return COALESCENT_HEIGHTS_REFERENCE_PRIOR;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            String fileName = xo.getStringAttribute(FileHelpers.FILE_NAME);

//            try {

//                File file = new File(fileName);
//                String parent = file.getParent();
//
//                if (!file.isAbsolute()) {
//                    parent = System.getProperty("user.dir");
//                }
//                file = new File(parent, fileName);
//                fileName = file.getAbsolutePath();
//
//                String parameterName = xo.getStringAttribute(PARAMETER_COLUMN);
//                int dimension = xo.getIntegerAttribute("dimension");
//
//                LogFileTraces traces = new LogFileTraces(fileName, file);
//                traces.loadTraces();
//                int maxState = traces.getMaxState();
//
//                // leaving the burnin attribute off will result in 10% being used
//                int burnin = xo.getAttribute("burnin", maxState / 10);
//                if (burnin < 0 || burnin >= maxState) {
//                    burnin = maxState / 10;
//                    System.out.println("WARNING: Burn-in larger than total number of states - using 10%");
//                }
//                traces.setBurnIn(burnin);
//
//                int[] traceIndexParameter = new int[dimension];
//                for (int i = 0; i < traceIndexParameter.length; i++) {
//					traceIndexParameter[i] = -1;
//				}
//
//                String[] columnNames = new String[dimension];
//				for (int i = 1; i <= columnNames.length; i++) {
//					columnNames[i-1] = parameterName + i;
//				}
//				System.out.println("Looking for the following columns:");
//				for (int i = 0; i < columnNames.length; i++) {
//					System.out.println("  " + columnNames[i]);
//				}
//				for (int i = 0; i < traces.getTraceCount(); i++) {
//					String traceName = traces.getTraceName(i);
//					for (int j = 0; j < columnNames.length; j++) {
//						if (traceName.trim().equals(columnNames[j])) {
//							traceIndexParameter[j] = i;
//							break;
//						}
//					}
//				}
//				System.out.println("Overview of traceIndexParameter:");
//				for (int i = 0; i < traceIndexParameter.length; i++) {
//					if (traceIndexParameter[i] == -1) {
//						throw new XMLParseException("Not all traces could be linked to the required columns.");
//					}
//					System.out.println("  traceIndexParameter[" + i + "]: " + traceIndexParameter[i] );
//				}
//
//				boolean[] flags = new boolean[dimension];
//				for (int i = 0; i < dimension; i++) {
//					flags[i] = true;
//				}
//
//				Double[][] parameterSamples = new Double[dimension][traces.getStateCount()];
//				for (int i = 0; i < dimension; i++) {
//					traces.getValues(traceIndexParameter[i]).toArray(parameterSamples[i]);
//					double initial = parameterSamples[i][0];
//					boolean tempFlag = false;
//					for (int j = 0; j < parameterSamples[i].length; j++) {
//						if (parameterSamples[i][j] != initial) {
//							tempFlag = true;
//							break;
//						}
//					}
//					flags[i] = tempFlag;
//				}

            final int dim = 8;

            boolean[] flags = new boolean[dim];
            double[] shapes = new double[dim];
            double[] scales = new double[dim];

            for (int i = 4; i < dim; ++i) flags[i] = true;


            /*for (int i = 0; i < flags.length; i++) {
                       if (flags[i]) {
                           double mean = 0.0;
                           for (int j = 0; j < parameterSamples[i].length; j++) {
                               mean += parameterSamples[i][j];
                           }
                           mean /= ((double)parameterSamples[i].length);
                           double variance = 0.0;
                           for (int j = 0; j < parameterSamples[i].length; j++) {
                               variance += Math.pow(parameterSamples[i][j] - mean, 2);
                           }
                           variance /= ((double)(parameterSamples[i].length-1));
                           scales[i] = variance/mean;
                           shapes[i] = mean/scales[i];
                           System.err.println(i + ": " + scales[i] + "   " + shapes[i]);
                       }
                   }*/

            for (int i = 0; i < flags.length; i++) {
                if (flags[i]) {
                    shapes[i] = 1.0;
                    scales[i] = 60.9194 / (double) (Binomial.choose2(9 - i));
                }
            }

//                System.err.println("4 choose 2 = " + Binomial.choose2(4));

            MultivariateGammaDistribution mvgd = new MultivariateGammaDistribution(shapes, scales, flags) {
                public double logPdf(double[] x) {
                    double rtn = super.logPdf(x);
//                        System.err.println("intervals = " + new Vector(x) + " " + rtn);
//                        System.exit(-1);
                    return rtn;
                }

                ;

            };

            MultivariateDistributionLikelihood likelihood = new MultivariateDistributionLikelihood(mvgd);
            for (int j = 0; j < xo.getChildCount(); j++) {
                if (xo.getChild(j) instanceof Statistic) {
                    likelihood.addData((Statistic) xo.getChild(j));
                } else {
                    throw new XMLParseException("illegal element in " + xo.getName() + " element");
                }
            }

            return likelihood;

//            } catch (FileNotFoundException fnfe) {
//                throw new XMLParseException("File '" + fileName + "' can not be opened for " + getParserName() + " element.");
//            } catch (java.io.IOException ioe) {
//                throw new XMLParseException(ioe.getMessage());
//            } catch (TraceException e) {
//                throw new XMLParseException(e.getMessage());
//            }

        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                AttributeRule.newStringRule("fileName"),
                AttributeRule.newStringRule("parameterColumn"),
                AttributeRule.newIntegerRule("dimension"),
                //AttributeRule.newIntegerRule("burnin"),
                new ElementRule(Statistic.class, 1, Integer.MAX_VALUE)
        };

        public String getParserDescription() {
            return "Calculates the coalescent height probabilities based on a sample of coalescent heights.";
        }

        public Class getReturnType() {
            return MultivariateDistributionLikelihood.class;
        }

    };

    /*private boolean isVariable(Double[] samples) {
        if (samples.length == 0) {
            throw new XMLParseException("Column length is zero.");
        }
        double initial = samples[0];
        int stop = samples.length;
        boolean flag = false;
        for (int i = 0; i < stop; i++) {
            if (samples[i] != initial) {
                flag = true;
                break;
            }
        }
        return flag;
    }*/

    /**
     * A special parser that reads a convenient short form of reference priors on trees.
     */
    public static XMLObjectParser CONDITIONAL_CLADE_REFERENCE_PRIOR_PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return CONDITIONAL_CLADE_REFERENCE_PRIOR;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            //Coalescent.TREEPRIOR = true;

            TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
            if (false) {
                String fileName = xo.getStringAttribute(FileHelpers.FILE_NAME);

                try {

                    File file = new File(fileName);
                    String name = file.getName();
                    String parent = file.getParent();

                    if (!file.isAbsolute()) {
                        parent = System.getProperty("user.dir");
                    }
                    file = new File(parent, fileName);
                    fileName = file.getAbsolutePath();

                    Reader reader = new FileReader(new File(parent, name));

                    // the burn-in is used as the number of trees discarded
                    int burnin = -1;
                    if (xo.hasAttribute(BURNIN)) {
                        // leaving the burnin attribute off will result in 10% being used
                        burnin = xo.getIntegerAttribute(BURNIN);
                    }

                    // the epsilon value which represents the number of occurrences for every not observed clade
                    double e = 1.0;
                    if (xo.hasAttribute(EPSILON)) {
                        // leaving the epsilon attribute off will result in 1.0 being used
                        e = xo.getDoubleAttribute(EPSILON);
                    }

                    TreeTrace trace = TreeTrace.loadTreeTrace(reader);

                    ConditionalCladeFrequency ccf = new ConditionalCladeFrequency(new TreeTrace[]{trace}, e, burnin, false);

                    ConditionalCladeProbability ccp = new ConditionalCladeProbability(ccf, treeModel);

                    return ccp;

                } catch (FileNotFoundException fnfe) {
                    throw new XMLParseException("File '" + fileName + "' can not be opened for " + getParserName() + " element.");
                } catch (java.io.IOException ioe) {
                    throw new XMLParseException(ioe.getMessage());
                } catch (ImportException ie) {
                    throw new XMLParseException(ie.getMessage());
                }
            }

            final int nTaxa = treeModel.getExternalNodeCount();
            double logPDF = 0.0;
            for (int i = nTaxa; i > 2; --i) {
                logPDF += Math.log(Binomial.choose2(i));
            }

//            System.err.println("logPDF = " + logPDF);
//            System.exit(-1);

            return new ConstantLikelihood(-logPDF);
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                AttributeRule.newStringRule("fileName"),
                AttributeRule.newIntegerRule(BURNIN),
                AttributeRule.newDoubleRule(EPSILON),
                new ElementRule(TreeModel.class)
        };

        public String getParserDescription() {
            return "Calculates the conditional clade probability of a tree based on a sample of tree space.";
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

            double mean = xo.getDoubleAttribute(MEAN);
            final double stdev = xo.getDoubleAttribute(STDEV);
            final double offset = xo.getAttribute(OFFSET, 0.0);
            final boolean meanInRealSpace = xo.getAttribute(MEAN_IN_REAL_SPACE, false);

            if (meanInRealSpace) {
                if (mean <= 0) {
                    throw new IllegalArgumentException("meanInRealSpace works only for a positive mean");
                }
                mean = Math.log(mean) - 0.5 * stdev * stdev;
            }

            final DistributionLikelihood likelihood = new DistributionLikelihood(new LogNormalDistribution(mean, stdev), offset);

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
            final double offset = xo.getDoubleAttribute(OFFSET);

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
                AttributeRule.newDoubleRule(OFFSET),
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

            DistributionLikelihood likelihood = new DistributionLikelihood(new BetaDistribution(shape, shapeB), offset);
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

}
