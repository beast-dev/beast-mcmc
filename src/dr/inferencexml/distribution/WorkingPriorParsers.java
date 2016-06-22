/*
 * WorkingPriorParsers.java
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

import dr.inference.distribution.DistributionLikelihood;
import dr.inference.distribution.MultivariateDistributionLikelihood;
import dr.inference.model.Likelihood;
import dr.inference.model.Statistic;
import dr.inference.trace.LogFileTraces;
import dr.inference.trace.TraceException;
import dr.math.distributions.*;
import dr.util.FileHelpers;
import dr.xml.*;

import java.io.File;
import java.io.FileNotFoundException;

/**
 * @author Guy Baele
 * @author Marc Suchard
 */
public class WorkingPriorParsers {
    public final static boolean DEBUG = true;

    public static final String NORMAL_REFERENCE_PRIOR = "normalReferencePrior";
    public static final String NORMAL_WORKING_PRIOR = "normalWorkingPrior";
    public static final String LOG_TRANSFORMED_NORMAL_REFERENCE_PRIOR = "logTransformedNormalReferencePrior";
    public static final String LOG_TRANSFORMED_NORMAL_WORKING_PRIOR = "logTransformedNormalWorkingPrior";
    public static final String LOGIT_TRANSFORMED_NORMAL_REFERENCE_PRIOR = "logitTransformedNormalReferencePrior";
    public static final String LOGIT_TRANSFORMED_NORMAL_WORKING_PRIOR = "logitTransformedNormalWorkingPrior";
    public static final String GAMMA_REFERENCE_PRIOR = "gammaReferencePrior";
    public static final String GAMMA_WORKING_PRIOR = "gammaWorkingPrior";
    public static final String PARAMETER_COLUMN = "parameterColumn";
    public static final String DIMENSION = "dimension";
    public static final String UPPERLIMIT = "upperLimit";

    /**
     * A special parser that reads a convenient short form of reference priors on parameters.
     */
    public static XMLObjectParser GAMMA_REFERENCE_PRIOR_PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return GAMMA_REFERENCE_PRIOR;
        }

        public String[] getParserNames() {
            return new String[]{getParserName(), GAMMA_WORKING_PRIOR};
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

                int dimension = 1;
                if (xo.hasAttribute(DIMENSION)) {
                    dimension = xo.getIntegerAttribute(DIMENSION);
                }
                if (dimension <= 0) {
                    throw new XMLParseException("Column '" + parameterName + "' has dimension smaller than 1.");
                }

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

                if (dimension == 1) {

                    int traceIndexParameter = -1;
                    for (int i = 0; i < traces.getTraceCount(); i++) {
                        String traceName = traces.getTraceName(i);
                        if (traceName.trim().equals(parameterName)) {
                            traceIndexParameter = i;
                        }
                    }

                    if (traceIndexParameter == -1) {
                        throw new XMLParseException("GammaKDEDistribution: Column '" + parameterName + "' can not be found for " + getParserName() + " element.");
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

                } else {

                    //dimension > 1
                    GammaKDEDistribution[] arrayKDE = new GammaKDEDistribution[dimension];

                    for (int i = 0; i < dimension; i++) {
                        //look for parameterName1, parameterName2, ... if necessary
                        String newParameterName = parameterName + (i+1);
                        int traceIndexParameter = -1;
                        for (int j = 0; j < traces.getTraceCount(); j++) {
                            String traceName = traces.getTraceName(j);
                            if (traceName.trim().equals(newParameterName)) {
                                traceIndexParameter = j;
                            }
                        }

                        if (traceIndexParameter == -1) {
                            throw new XMLParseException("GammaKDEDistribution: Column '" + newParameterName + "' can not be found for " + getParserName() + " element.");
                        }

                        Double[] parameterSamples = new Double[traces.getStateCount()];
                        traces.getValues(traceIndexParameter).toArray(parameterSamples);

                        arrayKDE[i] = new GammaKDEDistribution(parameterSamples);

                    }

                    MultivariateDistributionLikelihood likelihood = new MultivariateDistributionLikelihood(new MultivariateKDEDistribution(arrayKDE));

                    for (int j = 0; j < xo.getChildCount(); j++) {
                        if (xo.getChild(j) instanceof Statistic) {
                            if (DEBUG) {
                                System.out.println(((Statistic) xo.getChild(j)).toString());
                                System.out.println(((Statistic) xo.getChild(j)).getDimension());
                            }
                            likelihood.addData((Statistic) xo.getChild(j));
                        } else {
                            throw new XMLParseException("illegal element in " + xo.getName() + " element");
                        }
                    }

                    return likelihood;

                }

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

        public String[] getParserNames() {
            return new String[]{getParserName(), LOG_TRANSFORMED_NORMAL_WORKING_PRIOR};
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

                int dimension = 1;
                if (xo.hasAttribute(DIMENSION)) {
                    dimension = xo.getIntegerAttribute(DIMENSION);
                }
                if (dimension <= 0) {
                    throw new XMLParseException("Column '" + parameterName + "' has dimension smaller than 1.");
                }

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

                if (dimension == 1) {

                    int traceIndexParameter = -1;
                    for (int i = 0; i < traces.getTraceCount(); i++) {
                        String traceName = traces.getTraceName(i);
                        if (traceName.trim().equals(parameterName)) {
                            traceIndexParameter = i;
                        }
                    }

                    if (traceIndexParameter == -1) {
                        throw new XMLParseException("LogTransformedNormalKDEDistribution: Column '" + parameterName + "' can not be found for " + getParserName() + " element.");
                    }

                    Double[] parameterSamples = new Double[traces.getStateCount()];
                    traces.getValues(traceIndexParameter).toArray(parameterSamples);

                    DistributionLikelihood likelihood = new DistributionLikelihood(new LogTransformedNormalKDEDistribution(parameterSamples));
                    for (int j = 0; j < xo.getChildCount(); j++) {
                        if (xo.getChild(j) instanceof Statistic) {
                            if (DEBUG) {
                                System.out.println(((Statistic) xo.getChild(j)).toString());
                                System.out.println(((Statistic) xo.getChild(j)).getDimension());
                            }
                            likelihood.addData((Statistic) xo.getChild(j));
                        } else {
                            throw new XMLParseException("illegal element in " + xo.getName() + " element");
                        }
                    }

                    return likelihood;

                } else {

                    //dimension > 1
                    LogTransformedNormalKDEDistribution[] arrayKDE = new LogTransformedNormalKDEDistribution[dimension];

                    for (int i = 0; i < dimension; i++) {
                        //look for parameterName1, parameterName2, ... if necessary
                        String newParameterName = parameterName + (i+1);
                        int traceIndexParameter = -1;
                        for (int j = 0; j < traces.getTraceCount(); j++) {
                            String traceName = traces.getTraceName(j);
                            if (traceName.trim().equals(newParameterName)) {
                                traceIndexParameter = j;
                            }
                        }

                        if (traceIndexParameter == -1) {
                            throw new XMLParseException("LogTransformedNormalKDEDistribution: Column '" + newParameterName + "' can not be found for " + getParserName() + " element.");
                        }

                        Double[] parameterSamples = new Double[traces.getStateCount()];
                        traces.getValues(traceIndexParameter).toArray(parameterSamples);

                        arrayKDE[i] = new LogTransformedNormalKDEDistribution(parameterSamples);

                    }

                    MultivariateDistributionLikelihood likelihood = new MultivariateDistributionLikelihood(new MultivariateKDEDistribution(arrayKDE));

                    for (int j = 0; j < xo.getChildCount(); j++) {
                        if (xo.getChild(j) instanceof Statistic) {
                            if (DEBUG) {
                                System.out.println(((Statistic) xo.getChild(j)).toString());
                                System.out.println(((Statistic) xo.getChild(j)).getDimension());
                            }
                            likelihood.addData((Statistic) xo.getChild(j));
                        } else {
                            throw new XMLParseException("illegal element in " + xo.getName() + " element");
                        }
                    }

                    return likelihood;

                }

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
    public static XMLObjectParser LOGIT_TRANSFORMED_NORMAL_REFERENCE_PRIOR_PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return LOGIT_TRANSFORMED_NORMAL_REFERENCE_PRIOR;
        }

        public String[] getParserNames() {
            return new String[]{getParserName(), LOGIT_TRANSFORMED_NORMAL_WORKING_PRIOR};
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

                //keep using String and not an array of Strings, append integers later on
                String parameterName = xo.getStringAttribute(PARAMETER_COLUMN);

                int dimension = 1;
                if (xo.hasAttribute(DIMENSION)) {
                    dimension = xo.getIntegerAttribute(DIMENSION);
                }
                if (dimension <= 0) {
                    throw new XMLParseException("Column '" + parameterName + "' has dimension smaller than 1.");
                }

                double upperlimit = 1.0;
                if (xo.hasAttribute(UPPERLIMIT)) {
                    upperlimit = xo.getDoubleAttribute(UPPERLIMIT);
                }
                if (upperlimit <= 0.0) {
                    throw new XMLParseException("Positive upper bound expected for logit transformed normal KDE distribution.");
                }

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

                if (dimension == 1) {

                    int traceIndexParameter = -1;
                    for (int i = 0; i < traces.getTraceCount(); i++) {
                        String traceName = traces.getTraceName(i);
                        if (traceName.trim().equals(parameterName)) {
                            traceIndexParameter = i;
                        }
                    }

                    if (traceIndexParameter == -1) {
                        throw new XMLParseException("LogitTransformedNormalKDEDistribution: Column '" + parameterName + "' can not be found for " + getParserName() + " element.");
                    }

                    Double[] parameterSamples = new Double[traces.getStateCount()];
                    traces.getValues(traceIndexParameter).toArray(parameterSamples);

                    DistributionLikelihood likelihood = new DistributionLikelihood(new LogitTransformedNormalKDEDistribution(parameterSamples, upperlimit));
                    for (int j = 0; j < xo.getChildCount(); j++) {
                        if (xo.getChild(j) instanceof Statistic) {
                            if (DEBUG) {
                                System.out.println(((Statistic) xo.getChild(j)).toString());
                                System.out.println(((Statistic) xo.getChild(j)).getDimension());
                            }
                            likelihood.addData((Statistic) xo.getChild(j));
                        } else {
                            throw new XMLParseException("illegal element in " + xo.getName() + " element");
                        }
                    }

                    return likelihood;

                } else {

                    //dimension > 1
                    LogitTransformedNormalKDEDistribution[] arrayKDE = new LogitTransformedNormalKDEDistribution[dimension];

                    for (int i = 0; i < dimension; i++) {
                        //look for parameterName1, parameterName2, ... if necessary
                        String newParameterName = parameterName + (i+1);
                        int traceIndexParameter = -1;
                        for (int j = 0; j < traces.getTraceCount(); j++) {
                            String traceName = traces.getTraceName(j);
                            if (traceName.trim().equals(newParameterName)) {
                                traceIndexParameter = j;
                            }
                        }

                        if (traceIndexParameter == -1) {
                            throw new XMLParseException("LogitTransformedNormalKDEDistribution: Column '" + newParameterName + "' can not be found for " + getParserName() + " element.");
                        }

                        Double[] parameterSamples = new Double[traces.getStateCount()];
                        traces.getValues(traceIndexParameter).toArray(parameterSamples);

                        arrayKDE[i] = new LogitTransformedNormalKDEDistribution(parameterSamples, upperlimit);

                    }

                    MultivariateDistributionLikelihood likelihood = new MultivariateDistributionLikelihood(new MultivariateKDEDistribution(arrayKDE));

                    for (int j = 0; j < xo.getChildCount(); j++) {
                        if (xo.getChild(j) instanceof Statistic) {
                            if (DEBUG) {
                                System.out.println(((Statistic) xo.getChild(j)).toString());
                                System.out.println(((Statistic) xo.getChild(j)).getDimension());
                            }
                            likelihood.addData((Statistic) xo.getChild(j));
                        } else {
                            throw new XMLParseException("illegal element in " + xo.getName() + " element");
                        }
                    }

                    return likelihood;

                }

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
                //optional to provide a dimension attribute
                AttributeRule.newIntegerRule("dimension", true),
                //optional to provide an upperLimit attribute
                AttributeRule.newDoubleRule("upperLimit", true),
                new ElementRule(Statistic.class, 1, Integer.MAX_VALUE)
        };

        public String getParserDescription() {
            return "Calculates the reference prior probability of some data under logit transformed normal distribution.";
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

        public String[] getParserNames() {
            return new String[]{getParserName(), NORMAL_WORKING_PRIOR};
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

                int dimension = 1;
                if (xo.hasAttribute(DIMENSION)) {
                    dimension = xo.getIntegerAttribute(DIMENSION);
                }
                if (dimension <= 0) {
                    throw new XMLParseException("Column '" + parameterName + "' has dimension smaller than 1.");
                }

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

                if (dimension == 1) {

                    int traceIndexParameter = -1;
                    for (int i = 0; i < traces.getTraceCount(); i++) {
                        String traceName = traces.getTraceName(i);
                        if (traceName.trim().equals(parameterName)) {
                            traceIndexParameter = i;
                        }
                    }

                    if (traceIndexParameter == -1) {
                        throw new XMLParseException("NormalKDEDistribution: Column '" + parameterName + "' can not be found for " + getParserName() + " element.");
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

                } else {

                    //dimension > 1
                    NormalKDEDistribution[] arrayKDE = new NormalKDEDistribution[dimension];

                    for (int i = 0; i < dimension; i++) {
                        //look for parameterName1, parameterName2, ... if necessary
                        String newParameterName = parameterName + (i+1);
                        int traceIndexParameter = -1;
                        for (int j = 0; j < traces.getTraceCount(); j++) {
                            String traceName = traces.getTraceName(j);
                            if (traceName.trim().equals(newParameterName)) {
                                traceIndexParameter = j;
                            }
                        }

                        if (traceIndexParameter == -1) {
                            throw new XMLParseException("NormalKDEDistribution: Column '" + newParameterName + "' can not be found for " + getParserName() + " element.");
                        }

                        Double[] parameterSamples = new Double[traces.getStateCount()];
                        traces.getValues(traceIndexParameter).toArray(parameterSamples);

                        arrayKDE[i] = new NormalKDEDistribution(parameterSamples);

                    }

                    MultivariateDistributionLikelihood likelihood = new MultivariateDistributionLikelihood(new MultivariateKDEDistribution(arrayKDE));

                    for (int j = 0; j < xo.getChildCount(); j++) {
                        if (xo.getChild(j) instanceof Statistic) {
                            if (DEBUG) {
                                System.out.println(((Statistic) xo.getChild(j)).toString());
                                System.out.println(((Statistic) xo.getChild(j)).getDimension());
                            }
                            likelihood.addData((Statistic) xo.getChild(j));
                        } else {
                            throw new XMLParseException("illegal element in " + xo.getName() + " element");
                        }
                    }

                    return likelihood;

                }

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

}
