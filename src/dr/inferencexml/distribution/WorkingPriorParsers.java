/*
 * WorkingPriorParsers.java
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

import dr.inference.distribution.DistributionLikelihood;
import dr.inference.model.Likelihood;
import dr.inference.model.Statistic;
import dr.inference.trace.LogFileTraces;
import dr.inference.trace.TraceException;
import dr.math.distributions.GammaKDEDistribution;
import dr.math.distributions.LogTransformedNormalKDEDistribution;
import dr.math.distributions.LogitTransformedNormalKDEDistribution;
import dr.math.distributions.NormalKDEDistribution;
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
    public static final String LOG_TRANSFORMED_NORMAL_REFERENCE_PRIOR = "logTransformedNormalReferencePrior";
    public static final String LOGIT_TRANSFORMED_NORMAL_REFERENCE_PRIOR = "logitTransformedNormalReferencePrior";
    public static final String GAMMA_REFERENCE_PRIOR = "gammaReferencePrior";
    public static final String PARAMETER_COLUMN = "parameterColumn";

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
    public static XMLObjectParser LOGIT_TRANSFORMED_NORMAL_REFERENCE_PRIOR_PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return LOGIT_TRANSFORMED_NORMAL_REFERENCE_PRIOR;
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

                DistributionLikelihood likelihood = new DistributionLikelihood(new LogitTransformedNormalKDEDistribution(parameterSamples));
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

}
