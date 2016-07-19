/*
 * TreeWorkingPriorParsers.java
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

package dr.evomodelxml;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.Reader;

import jebl.math.Binomial;

import org.apache.commons.math.MathException;
import org.apache.commons.math.analysis.interpolation.LoessInterpolator;

import dr.evolution.io.Importer.ImportException;
import dr.evolution.io.TreeTrace;
import dr.evomodel.coalescent.CoalescentConstantLikelihood;
import dr.evomodel.coalescent.ExponentialProductLikelihood;
import dr.evomodel.coalescent.ExponentialProductPosteriorMeansLikelihood;
import dr.evomodel.coalescent.ExponentialProductSufficientStatisticsLikelihood;
import dr.evomodel.coalescent.GammaProductLikelihood;
import dr.evomodel.tree.ConditionalCladeFrequency;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.tree.ConditionalCladeProbability;
import dr.inference.distribution.MultivariateDistributionLikelihood;
import dr.inference.model.ConstantLikelihood;
import dr.inference.model.Likelihood;
import dr.inference.model.Statistic;
import dr.inference.trace.LogFileTraces;
import dr.inference.trace.TraceException;
import dr.math.distributions.MultivariateGammaDistribution;
import dr.math.functionEval.GammaFunction;
import dr.util.FileHelpers;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.AttributeRule;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLObjectParser;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

/**
 * @author Guy Baele
 * @author Marc A. Suchard
 */
public class TreeWorkingPriorParsers {

    public final static boolean DEBUG = false;

    public static final String CONSTANT_TREE_TOPOLOGY_PRIOR = "constantTreeTopologyPrior";
    public static final String CONTEMPORANEOUS_COALESCENT_CONSTANT = "contemporaneousCoalescentConstantLikelihood";
    public static final String COALESCENT_CONSTANT_LIKELIHOOD = "coalescentConstantLikelihood";
    public static final String CONDITIONAL_CLADE_REFERENCE_PRIOR = "conditionalCladeProbability";
    public static final String COALESCENT_HEIGHTS_REFERENCE_PRIOR = "coalescentHeightsReferencePrior";
    public static final String PRODUCT_OF_EXPONENTIALS = "productOfExponentials";
    public static final String PRODUCT_OF_EXPONENTIALS_POSTERIOR_MEANS = "productOfExponentialsPosteriorMeans";
    public static final String PRODUCT_OF_EXPONENTIALS_POSTERIOR_MEANS_LOESS = "productOfExponentialsPosteriorMeansLoess";
    public static final String PRODUCT_OF_EXPONENTIALS_SUFFICIENT_STATISTICS = "productOfExponentialsSufficientStatistics";
    public static final String PRODUCT_OF_GAMMAS = "productOfGammas";
    public static final String CONSTANT_DECREASED_VARIANCE_PRIOR = "constantDecreasedVariancePrior";
    public final static String BURNIN = "burnin";
    public final static String EPSILON = "epsilon";
    public static final String PARAMETER_COLUMN = "parameterColumn";

    public static final double epsilon = 0.1;

    /**
     * 
     */
    public static XMLObjectParser COALESCENT_CONSTANT_LIKELIHOOD_PARSER = new AbstractXMLObjectParser () {

        public String getParserName() {
            return COALESCENT_CONSTANT_LIKELIHOOD;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            //only use as test class, will require separate likelihood function
            /*Tree tree = (Tree) xo.getChild(Tree.class);
    		TreeIntervals intervals = new TreeIntervals(tree);
    		final int nIntervals = intervals.getIntervalCount();
    		System.err.println("Interval count: " + nIntervals);
    		double logPDF = 0.0;
            for (int i = 0; i < nIntervals; i++) {
            	System.err.println("Lineage count " + i + ": " + intervals.getLineageCount(i));
            	if (intervals.getLineageCount(i) > 2) {
            		logPDF += Math.log(Binomial.choose2(intervals.getLineageCount(i)));
            	}
            }
            System.err.println("logPDF = " + logPDF);
            return new ConstantLikelihood(-logPDF);*/

            TreeModel tree = (TreeModel) xo.getChild(TreeModel.class);
            //System.err.println(tree);
            return new CoalescentConstantLikelihood(tree);

        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                new ElementRule(TreeModel.class)
        };

        public String getParserDescription() {
            return "Calculates the number of possible combinations of coalescent events.";
        }

        public Class getReturnType() {
            return Likelihood.class;
        }

    };

    /**
     * 
     */
    public static XMLObjectParser CONTEMPORANEOUS_COALESCENT_CONSTANT_PARSER = new AbstractXMLObjectParser () {

        public String getParserName() {
            return CONTEMPORANEOUS_COALESCENT_CONSTANT;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
            final int nTaxa = treeModel.getExternalNodeCount();
            double logPDF = 0.0;
            for (int i = nTaxa; i > 2; --i) {
                logPDF += Math.log(Binomial.choose2(i));
            }
            return new ConstantLikelihood(-logPDF);

        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                new ElementRule(TreeModel.class)
        };

        public String getParserDescription() {
            return "Calculates the number of possible combinations of coalescent events.";
        }

        public Class getReturnType() {
            return Likelihood.class;
        }

    };

    /**
     * Parser that aims to mimic the behaviour of the GMRF Skyride for fixed log population sizes.
     */
    public static XMLObjectParser PRODUCT_OF_EXPONENTIALS_PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return PRODUCT_OF_EXPONENTIALS;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            double logPopSize = xo.getDoubleAttribute("logPopSize");
            TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);

            return new ExponentialProductLikelihood(treeModel, logPopSize);

        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                AttributeRule.newDoubleRule("logPopSize"),
                new ElementRule(TreeModel.class)
        };

        public String getParserDescription() {
            return "Calculates a product of exponentials based on a (set of) fixed log population sizes.";
        }

        public Class getReturnType() {
            return Likelihood.class;
        }

    };

    /**
     * Parser that aims to use the posterior means of the time between coalescent events.
     */
    public static XMLObjectParser CONSTANT_DECREASED_VARIANCE_PRIOR_PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return CONSTANT_DECREASED_VARIANCE_PRIOR;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
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
                int dimension = xo.getIntegerAttribute("dimension");

                LogFileTraces traces = new LogFileTraces(fileName, file);
                traces.loadTraces();
                long maxState = traces.getMaxState();

                // leaving the burnin attribute off will result in 10% being used
                long burnin = xo.getAttribute("burnin", maxState / 10);
                if (burnin < 0 || burnin >= maxState) {
                    burnin = maxState / 10;
                    System.out.println("WARNING: Burn-in larger than total number of states - using 10%");
                }
                traces.setBurnIn(burnin);

                int traceIndexParameter = -1;

                if (DEBUG) {
                    System.err.println("Looking for the following column:" + parameterName);
                }
                for (int i = 0; i < traces.getTraceCount(); i++) {
                    String traceName = traces.getTraceName(i);
                    if (traceName.trim().equals(parameterName)) {
                        traceIndexParameter = i;
                        break;
                    }
                }
                if (DEBUG) {
                    System.err.println("Overview of traceIndexParameter:");
                }
                if (traceIndexParameter == -1) {
                    throw new XMLParseException("Not all traces could be linked to the required columns.");
                } else if (DEBUG) {
                    System.out.println("  traceIndexParameter: " + traceIndexParameter);
                }

                Double[] parameterSamples = new Double[traces.getStateCount()];
                traces.getValues(traceIndexParameter).toArray(parameterSamples);

                //not necessary to work with flags, assume that we are using coalescentEventsStatistic
                double posteriorMean = 0.0;
                for (int i = 0; i < parameterSamples.length; i++) {
                    posteriorMean += parameterSamples[i];
                }
                posteriorMean /= ((double)parameterSamples.length);
                if (DEBUG) {
                    System.err.println("Variable column -> " + posteriorMean);
                }
                //posteriorMean = Math.log(posteriorMean);
                //System.err.println("Log transformed: " + posteriorMean);

                //return new ExponentialProductPosteriorMeansLikelihood(treeModel, posteriorMeans);
                return new GammaProductLikelihood(treeModel, posteriorMean, dimension);

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
                AttributeRule.newIntegerRule("dimension"),
                new ElementRule(TreeModel.class)
        };

        public String getParserDescription() {
            return "Calculates a product of exponentials based on a set of posterior sample means.";
        }

        public Class getReturnType() {
            return Likelihood.class;
        }

    }; 

    /**
     * Parser that aims to use the posterior means of the time between coalescent events.
     */
    public static XMLObjectParser PRODUCT_OF_EXPONENTIALS_POSTERIOR_MEANS_PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return PRODUCT_OF_EXPONENTIALS_POSTERIOR_MEANS;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
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
                int dimension = xo.getIntegerAttribute("dimension");

                LogFileTraces traces = new LogFileTraces(fileName, file);
                traces.loadTraces();
                long maxState = traces.getMaxState();

                // leaving the burnin attribute off will result in 10% being used
                long burnin = xo.getAttribute("burnin", maxState / 10);
                if (burnin < 0 || burnin >= maxState) {
                    burnin = maxState / 10;
                    System.out.println("WARNING: Burn-in larger than total number of states - using 10%");
                }
                traces.setBurnIn(burnin);

                int[] traceIndexParameter = new int[dimension];
                for (int i = 0; i < traceIndexParameter.length; i++) {
                    traceIndexParameter[i] = -1;
                }

                String[] columnNames = new String[dimension];
                for (int i = 1; i <= columnNames.length; i++) {
                    columnNames[i-1] = parameterName + i;
                }
                if (DEBUG) {
                    System.err.println("Looking for the following columns:");
                    for (int i = 0; i < columnNames.length; i++) {
                        System.err.println("  " + columnNames[i]);
                    }
                }
                for (int i = 0; i < traces.getTraceCount(); i++) {
                    String traceName = traces.getTraceName(i);
                    for (int j = 0; j < columnNames.length; j++) {
                        if (traceName.trim().equals(columnNames[j])) {
                            traceIndexParameter[j] = i;
                            break;
                        }
                    }
                }
                if (DEBUG) {
                    System.err.println("Overview of traceIndexParameter:");
                }
                for (int i = 0; i < traceIndexParameter.length; i++) {
                    if (traceIndexParameter[i] == -1) {
                        throw new XMLParseException("Not all traces could be linked to the required columns, problem with trace index: " + i + "; traceIndexParameter.length = " + traceIndexParameter.length);
                    } else if (DEBUG) {
                        System.err.println("  traceIndexParameter[" + i + "]: " + traceIndexParameter[i]);
                    }
                }

                Double[][] parameterSamples = new Double[dimension][traces.getStateCount()];
                for (int i = 0; i < dimension; i++) {
                    traces.getValues(traceIndexParameter[i]).toArray(parameterSamples[i]);
                }

                //not necessary to work with flags, assume that we are using coalescentEventsStatistic
                double[] posteriorMeans = new double[dimension];

                for (int i = 0; i < posteriorMeans.length; i++) {
                    //variable column
                    double mean = 0.0;
                    for (int j = 0; j < parameterSamples[i].length; j++) {
                        mean += parameterSamples[i][j];
                    }
                    mean /= ((double)parameterSamples[i].length);

                    //mean = 1.0/mean;

                    if (DEBUG) {
                        System.err.println("Variable column: " + i + " -> " + mean);
                    }
                    posteriorMeans[i] = Math.log(mean);
                }

                return new ExponentialProductPosteriorMeansLikelihood(treeModel, posteriorMeans);

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
                AttributeRule.newIntegerRule("dimension"),
                new ElementRule(TreeModel.class)
        };

        public String getParserDescription() {
            return "Calculates a product of exponentials based on a set of posterior sample means.";
        }

        public Class getReturnType() {
            return Likelihood.class;
        }

    }; 

    /**
     * Parser that aims to use the posterior means of the time between coalescent events.
     */
    public static XMLObjectParser PRODUCT_OF_EXPONENTIALS_POSTERIOR_MEANS_LOESS_PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return PRODUCT_OF_EXPONENTIALS_POSTERIOR_MEANS_LOESS;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
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
                int dimension = xo.getIntegerAttribute("dimension");

                LogFileTraces traces = new LogFileTraces(fileName, file);
                traces.loadTraces();
                long maxState = traces.getMaxState();

                // leaving the burnin attribute off will result in 10% being used
                long burnin = xo.getAttribute("burnin", maxState / 10);
                if (burnin < 0 || burnin >= maxState) {
                    burnin = maxState / 10;
                    System.out.println("WARNING: Burn-in larger than total number of states - using 10%");
                    System.out.println("Burnin: " + burnin);
                }
                traces.setBurnIn(burnin);

                int[] traceIndexParameter = new int[dimension];
                for (int i = 0; i < traceIndexParameter.length; i++) {
                    traceIndexParameter[i] = -1;
                }

                String[] columnNames = new String[dimension];
                for (int i = 1; i <= columnNames.length; i++) {
                    columnNames[i-1] = parameterName + i;
                }
                if (DEBUG) {
                    System.err.println("Looking for the following columns:");
                    for (int i = 0; i < columnNames.length; i++) {
                        System.err.println("  " + columnNames[i]);
                    }
                }
                for (int i = 0; i < traces.getTraceCount(); i++) {
                    String traceName = traces.getTraceName(i);
                    for (int j = 0; j < columnNames.length; j++) {
                        if (traceName.trim().equals(columnNames[j])) {
                            traceIndexParameter[j] = i;
                            break;
                        }
                    }
                }
                if (DEBUG) {
                    System.err.println("Overview of traceIndexParameter:");
                }
                for (int i = 0; i < traceIndexParameter.length; i++) {
                    if (traceIndexParameter[i] == -1) {
                        throw new XMLParseException("Not all traces could be linked to the required columns, problem with trace index: " + i + "; traceIndexParameter.length = " + traceIndexParameter.length);
                    }
                    if (DEBUG) {
                        System.err.println("  traceIndexParameter[" + i + "]: " + traceIndexParameter[i]);
                    }
                }

                Double[][] parameterSamples = new Double[dimension][traces.getStateCount()];
                for (int i = 0; i < dimension; i++) {
                    traces.getValues(traceIndexParameter[i]).toArray(parameterSamples[i]);
                }

                //not necessary to work with flags, assume that we are using coalescentEventsStatistic
                double[] posteriorMeans = new double[dimension];

                for (int i = 0; i < posteriorMeans.length; i++) {
                    //variable column
                    double mean = 0.0;
                    for (int j = 0; j < parameterSamples[i].length; j++) {
                        mean += parameterSamples[i][j];
                    }
                    mean /= ((double)parameterSamples[i].length);

                    //mean = 1.0/mean;

                    if (DEBUG) {
                        System.err.println("Variable column: " + i + " -> " + mean);
                    }
                    posteriorMeans[i] = Math.log(mean);
                }

                try {
                    //Print log posterior means
                    if (DEBUG) {
                        System.err.println("Log Posterior Means:");
                        for (int i = 0; i < posteriorMeans.length; i++) {
                            System.err.println(posteriorMeans[i]);
                        }
                    }
                    //Call Loess interpolator here
                    LoessInterpolator loess = new LoessInterpolator(1.0, 2);
                    double[] xvalues = new double[posteriorMeans.length];
                    for (int i = 0; i < posteriorMeans.length; i++) {
                        xvalues[i] = i;
                    }

                    double[] loessOutput = loess.smooth(xvalues, posteriorMeans);

                    if (DEBUG) {
                        System.err.println("Loess output:");
                        for (int i = 0; i < loessOutput.length; i++) {
                            System.err.println(loessOutput[i]);
                        }
                    }

                    posteriorMeans = loessOutput;

                } catch (MathException e) {
                    // TODO Auto-generated catch block
                    e.printStackTrace();
                }

                return new ExponentialProductPosteriorMeansLikelihood(treeModel, posteriorMeans);

                /*for (int i = 0; i < posteriorMeans.length; i++) {
					posteriorMeans[i] = Math.exp(posteriorMeans[i]);
				}
				return new GammaIntervalProductLikelihood(treeModel, 0.0, posteriorMeans, posteriorMeans);*/

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
                AttributeRule.newIntegerRule("dimension"),
                new ElementRule(TreeModel.class)
        };

        public String getParserDescription() {
            return "Calculates a product of exponentials based on a set of posterior sample means, undergoing loess smoothing.";
        }

        public Class getReturnType() {
            return Likelihood.class;
        }

    }; 

    /**
     * Parser that aims to use the posterior means of the time between coalescent events, using sufficient statistics.
     */
    public static XMLObjectParser PRODUCT_OF_EXPONENTIALS_SUFFICIENT_STATISTICS_PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return PRODUCT_OF_EXPONENTIALS_SUFFICIENT_STATISTICS;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
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
                int dimension = xo.getIntegerAttribute("dimension");

                LogFileTraces traces = new LogFileTraces(fileName, file);
                traces.loadTraces();
                long maxState = traces.getMaxState();

                // leaving the burnin attribute off will result in 10% being used
                long burnin = xo.getAttribute("burnin", maxState / 10);
                if (burnin < 0 || burnin >= maxState) {
                    burnin = maxState / 10;
                    System.out.println("WARNING: Burn-in larger than total number of states - using 10%");
                }
                traces.setBurnIn(burnin);

                int[] traceIndexParameter = new int[dimension];
                for (int i = 0; i < traceIndexParameter.length; i++) {
                    traceIndexParameter[i] = -1;
                }

                String[] columnNames = new String[dimension];
                for (int i = 1; i <= columnNames.length; i++) {
                    columnNames[i-1] = parameterName + i;
                }
                if (DEBUG) {
                    System.err.println("Looking for the following columns:");
                    for (int i = 0; i < columnNames.length; i++) {
                        System.err.println("  " + columnNames[i]);
                    }
                }
                for (int i = 0; i < traces.getTraceCount(); i++) {
                    String traceName = traces.getTraceName(i);
                    for (int j = 0; j < columnNames.length; j++) {
                        if (traceName.trim().equals(columnNames[j])) {
                            traceIndexParameter[j] = i;
                            break;
                        }
                    }
                }
                if (DEBUG) {
                    System.err.println("Overview of traceIndexParameter:");
                }
                for (int i = 0; i < traceIndexParameter.length; i++) {
                    if (traceIndexParameter[i] == -1) {
                        throw new XMLParseException("Not all traces could be linked to the required columns.");
                    }
                    if (DEBUG) {
                        System.err.println("  traceIndexParameter[" + i + "]: " + traceIndexParameter[i]);
                    }
                }

                Double[][] parameterSamples = new Double[dimension][traces.getStateCount()];
                for (int i = 0; i < dimension; i++) {
                    traces.getValues(traceIndexParameter[i]).toArray(parameterSamples[i]);
                }

                //not necessary to work with flags, assume that we are using coalescentEventsStatistic
                double[] posteriorMeans = new double[dimension];

                for (int i = 0; i < posteriorMeans.length; i++) {
                    //variable column
                    double mean = 0.0;
                    for (int j = 0; j < parameterSamples[i].length; j++) {
                        mean += parameterSamples[i][j];
                    }
                    mean /= ((double)parameterSamples[i].length);
                    if (DEBUG) {
                        System.err.println("Variable column: " + i + " -> " + mean);
                    }
                    posteriorMeans[i] = Math.log(mean);
                }

                return new ExponentialProductSufficientStatisticsLikelihood(treeModel, posteriorMeans);

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
                AttributeRule.newIntegerRule("dimension"),
                new ElementRule(TreeModel.class)
        };

        public String getParserDescription() {
            return "Calculates a product of exponentials based on a set of posterior sample means, using sufficient statistics.";
        }

        public Class getReturnType() {
            return Likelihood.class;
        }

    }; 

    /**
     * Parser that forms a product of gamma distributions to decrease model testing variance 
     */
    public static XMLObjectParser PRODUCT_OF_GAMMAS_PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return PRODUCT_OF_GAMMAS;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
            double popSize = xo.getDoubleAttribute("popSize");
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
                int dimension = xo.getIntegerAttribute("dimension");

                LogFileTraces traces = new LogFileTraces(fileName, file);
                traces.loadTraces();
                long maxState = traces.getMaxState();

                // leaving the burnin attribute off will result in 10% being used
                long burnin = xo.getAttribute("burnin", maxState / 10);
                if (burnin < 0 || burnin >= maxState) {
                    burnin = maxState / 10;
                    System.out.println("WARNING: Burn-in larger than total number of states - using 10%");
                }
                traces.setBurnIn(burnin);

                int[] traceIndexParameter = new int[dimension];
                for (int i = 0; i < traceIndexParameter.length; i++) {
                    traceIndexParameter[i] = -1;
                }

                String[] columnNames = new String[dimension];
                for (int i = 1; i <= columnNames.length; i++) {
                    columnNames[i-1] = parameterName + i;
                }
                if (DEBUG) {
                    System.err.println("Looking for the following columns:");
                    for (int i = 0; i < columnNames.length; i++) {
                        System.err.println("  " + columnNames[i]);
                    }
                }
                for (int i = 0; i < traces.getTraceCount(); i++) {
                    String traceName = traces.getTraceName(i);
                    for (int j = 0; j < columnNames.length; j++) {
                        if (traceName.trim().equals(columnNames[j])) {
                            traceIndexParameter[j] = i;
                            break;
                        }
                    }
                }
                if (DEBUG) {
                    System.err.println("Overview of traceIndexParameter:");
                }
                for (int i = 0; i < traceIndexParameter.length; i++) {
                    if (traceIndexParameter[i] == -1) {
                        throw new XMLParseException("Not all traces could be linked to the required columns.");
                    }
                    if (DEBUG) {
                        System.err.println("  traceIndexParameter[" + i + "]: " + traceIndexParameter[i]);
                    }
                }

                boolean[] flags = new boolean[dimension];
                for (int i = 0; i < dimension; i++) {
                    flags[i] = true;
                }

                //this is code to identify constant columns
                Double[][] parameterSamples = new Double[dimension][traces.getStateCount()];
                for (int i = 0; i < dimension; i++) {
                    traces.getValues(traceIndexParameter[i]).toArray(parameterSamples[i]);
                    double initial = parameterSamples[i][0];
                    boolean tempFlag = false;
                    for (int j = 0; j < parameterSamples[i].length; j++) {
                        if (parameterSamples[i][j] != initial) {
                            tempFlag = true;
                            break;
                        }
                    }
                    flags[i] = tempFlag;
                }

                //double[] shapes = new double[dimension];
                //double[] scales = new double[dimension];
                double[] means = new double[dimension];
                double[] variances = new double[dimension];

                for (int i = 0; i < flags.length; i++) {
                    if (flags[i]) {
                        //variable column
                        double mean = 0.0;
                        for (int j = 0; j < parameterSamples[i].length; j++) {
                            mean += parameterSamples[i][j];
                        }
                        mean /= ((double)parameterSamples[i].length);
                        means[i] = mean;
                        double variance = 0.0;
                        for (int j = 0; j < parameterSamples[i].length; j++) {
                            variance += Math.pow(parameterSamples[i][j] - mean, 2);
                        }
                        variance /= ((double)(parameterSamples[i].length-1));
                        variances[i] = variance;
                        //scales[i] = variance/mean;
                        //shapes[i] = mean/scales[i];
                        //System.err.println("Variable column: " + i + " -> " + shapes[i] + "   " + scales[i]);
                        if (DEBUG) {
                            System.err.println("Variable column: " + i + " -> " + means[i] + "   " + variances[i]);
                        }
                    } else {
                        //constant column
                        double mean = 0.0;
                        for (int j = 0; j < parameterSamples[i].length; j++) {
                            mean += parameterSamples[i][j];
                        }
                        mean /= ((double)parameterSamples[i].length);
                        means[i] = mean;
                        //variance will be 0.0, so add epsilon
                        variances[i] = 0.0;
                        variances[i] = epsilon;
                        //System.err.println("mean = " + mean + "   variance = " + variance);
                        //scales[i] = variance/mean;
                        //shapes[i] = mean/scales[i];
                        if (DEBUG) {
                            System.err.println("Constant column: " + i + " -> " + means[i] + "   " + variances[i]);
                        }
                    } 
                }

                //only provide the actual intervals, from the root downwards
                /*int newLength = 0;
				for (int i = 0; i < flags.length; i++) {
					if (means[i] == 0.0 && variances[i] == epsilon) {
						//do nothing at the moment
					} else {
						newLength++;
					}
				}
				System.err.println("Flags length = " + flags.length);
				System.err.println("New length = " + newLength);
				double[] newMeans = new double[newLength];
				double[] newVariances = new double[newLength];
				int counter = 0;
				for (int i = 0; i < flags.length; i++) {
					if (means[i] == 0.0 && variances[i] == epsilon) {
						//do nothing at the moment
					} else {
						newMeans[counter] = means[i];
						newVariances[counter] = variances[i];
						counter++;
					}
				}*/

                //return new GammaProductLikelihood(treeModel, popSize, newMeans, newVariances);
                //return new GammaIntervalProductLikelihood(treeModel, popSize, means, variances);
                //return new GammaIntervalProductLikelihood(treeModel, popSize, newMeans, newVariances);


                return new GammaProductLikelihood(treeModel, popSize, means, variances);


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
                //AttributeRule.newDoubleRule("popSize"),
                AttributeRule.newStringRule("fileName"),
                AttributeRule.newStringRule("parameterColumn"),
                AttributeRule.newIntegerRule("dimension"),
                AttributeRule.newDoubleRule("popSize"),
                //new ElementRule(Parameter.class),
                new ElementRule(TreeModel.class)
        };

        public String getParserDescription() {
            return "Calculates a product of gamma distributions.";
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

            try {

                File file = new File(fileName);
                String parent = file.getParent();

                if (!file.isAbsolute()) {
                    parent = System.getProperty("user.dir");
                }
                file = new File(parent, fileName);
                fileName = file.getAbsolutePath();

                String parameterName = xo.getStringAttribute(PARAMETER_COLUMN);
                int dimension = xo.getIntegerAttribute("dimension");

                LogFileTraces traces = new LogFileTraces(fileName, file);
                traces.loadTraces();
                long maxState = traces.getMaxState();

                // leaving the burnin attribute off will result in 10% being used
                long burnin = xo.getAttribute("burnin", maxState / 10);
                if (burnin < 0 || burnin >= maxState) {
                    burnin = maxState / 10;
                    System.out.println("WARNING: Burn-in larger than total number of states - using 10%");
                }
                traces.setBurnIn(burnin);

                int[] traceIndexParameter = new int[dimension];
                for (int i = 0; i < traceIndexParameter.length; i++) {
                    traceIndexParameter[i] = -1;
                }

                String[] columnNames = new String[dimension];
                for (int i = 1; i <= columnNames.length; i++) {
                    columnNames[i-1] = parameterName + i;
                }
                if (DEBUG) {
                    System.err.println("Looking for the following columns:");
                    for (int i = 0; i < columnNames.length; i++) {
                        System.err.println("  " + columnNames[i]);
                    }
                }
                for (int i = 0; i < traces.getTraceCount(); i++) {
                    String traceName = traces.getTraceName(i);
                    for (int j = 0; j < columnNames.length; j++) {
                        if (traceName.trim().equals(columnNames[j])) {
                            traceIndexParameter[j] = i;
                            break;
                        }
                    }
                }
                if (DEBUG) {
                    System.err.println("Overview of traceIndexParameter:");
                }
                for (int i = 0; i < traceIndexParameter.length; i++) {
                    if (traceIndexParameter[i] == -1) {
                        throw new XMLParseException("Not all traces could be linked to the required columns.");
                    }
                    if (DEBUG) {
                        System.err.println("  traceIndexParameter[" + i + "]: " + traceIndexParameter[i]);
                    }
                }

                boolean[] flags = new boolean[dimension];
                for (int i = 0; i < dimension; i++) {
                    flags[i] = true;
                }

                Double[][] parameterSamples = new Double[dimension][traces.getStateCount()];
                for (int i = 0; i < dimension; i++) {
                    traces.getValues(traceIndexParameter[i]).toArray(parameterSamples[i]);
                    double initial = parameterSamples[i][0];
                    boolean tempFlag = false;
                    for (int j = 0; j < parameterSamples[i].length; j++) {
                        if (parameterSamples[i][j] != initial) {
                            tempFlag = true;
                            break;
                        }
                    }
                    flags[i] = tempFlag;
                }

                double[] shapes = new double[dimension];
                double[] scales = new double[dimension];

                for (int i = 0; i < flags.length; i++) {
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
                        if (DEBUG) {
                            System.err.println("Variable column: " + i + " -> " + shapes[i] + "   " + scales[i]);
                        }
                    } /*else {
                    	   double mean = 0.0;
                           for (int j = 0; j < parameterSamples[i].length; j++) {
                               mean += parameterSamples[i][j];
                           }
                           mean /= ((double)parameterSamples[i].length);
                           //variance will be 0.0, so add epsilon
                           double variance = 0.000001;
                           System.err.println("mean = " + mean + "   variance = " + variance);
                           scales[i] = variance/mean;
                           shapes[i] = mean/scales[i];
                    	   System.err.println("Constant column: " + i + " -> " + shapes[i] + "   " + scales[i]);
                       } */
                }

                if (DEBUG) {
                    System.err.println("Columns to be evaluated:");
                    for (int i = 0; i < flags.length; i++) {
                        if (flags[i]) {
                            System.err.println("Column " + i);
                        }
                    }
                }

                MultivariateGammaDistribution mvgd = new MultivariateGammaDistribution(shapes, scales, flags);

                MultivariateDistributionLikelihood likelihood = new MultivariateDistributionLikelihood(mvgd);
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

    /**
     * A simple parser to provide a constant prior on the tree topology
     */

    public static XMLObjectParser CONSTANT_TREE_TOPOLOGY_PRIOR_PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return CONSTANT_TREE_TOPOLOGY_PRIOR;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
            final int nTaxa = treeModel.getExternalNodeCount();
            double combinations = GammaFunction.factorial(2*nTaxa-3)/(GammaFunction.factorial(nTaxa-2)*Math.pow(2.0,nTaxa-2));
            double logPDF = Math.log(1.0/combinations);
            return new ConstantLikelihood(logPDF);

        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                new ElementRule(TreeModel.class)
        };

        public String getParserDescription() {
            return "Calculates the constant tree topology prior, i.e. 1 over the total number of rooted bifurcating trees.";
        }

        public Class getReturnType() {
            return Likelihood.class;
        }

    };

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

            /*final int nTaxa = treeModel.getExternalNodeCount();
            double logPDF = 0.0;
            for (int i = nTaxa; i > 2; --i) {
                logPDF += Math.log(Binomial.choose2(i));
            }

            return new ConstantLikelihood(-logPDF);*/
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

}
