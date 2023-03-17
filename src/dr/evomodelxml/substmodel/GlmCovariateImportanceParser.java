/*
 * GlmCovariateImportanceParser.java
 *
 * Copyright (c) 2002-2023 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodelxml.substmodel;

import dr.evomodel.coalescent.OldGMRFSkyrideLikelihood;
import dr.evomodel.substmodel.GlmCovariateImportance;
import dr.evomodel.substmodel.OldGLMSubstitutionModel;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.inference.distribution.DistributionLikelihood;
import dr.inference.distribution.MultivariateDistributionLikelihood;
import dr.inference.distribution.shrinkage.OldBayesianBridgeLikelihood;
import dr.inference.loggers.Loggable;
import dr.inference.model.Statistic;
import dr.inference.trace.LogFileTraces;
import dr.inference.trace.TraceException;
import dr.math.distributions.GammaKDEDistribution;
import dr.math.distributions.MultivariateKDEDistribution;
import dr.util.FileHelpers;
import dr.xml.*;

import java.io.File;
import java.io.FileNotFoundException;

/**
 * @author Marc Suchard
 * @author Philippe Lemey
 */
public class GlmCovariateImportanceParser extends AbstractXMLObjectParser {
    public final static boolean DEBUG = true;

    private static final String PARSER_NAME = "glmCovariateImportance";
    private static final String PARAMETER_COLUMN = "parameterColumn";
    private static final String DIMENSION = "dimension";

    public String getParserName() {
        return PARSER_NAME;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        TreeDataLikelihood likelihood = (TreeDataLikelihood) xo.getChild(TreeDataLikelihood.class);
        OldGLMSubstitutionModel substitutionModel = (OldGLMSubstitutionModel) xo.getChild(OldBayesianBridgeLikelihood.class);

        return new GlmCovariateImportance(likelihood, substitutionModel);
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newStringRule("fileName", true),
            AttributeRule.newStringRule("parameterColumn", true),
            AttributeRule.newIntegerRule("burnIn", true),
            new ElementRule(TreeDataLikelihood.class),
            new ElementRule(OldGLMSubstitutionModel.class),
    };

    public String getParserDescription() {
        return "Calculates model deviance for each fixed effect in a phylogeographic GLM";
    }

    public Class getReturnType() {
        return Loggable.class;
    }

    private Object readLog(XMLObject xo) throws XMLParseException {

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
                    String newParameterName = parameterName + (i + 1);
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

        } catch (
                FileNotFoundException fnfe) {
            throw new XMLParseException("File '" + fileName + "' can not be opened for " + getParserName() + " element.");
        } catch (
                java.io.IOException ioe) {
            throw new XMLParseException(ioe.getMessage());
        } catch (
                TraceException e) {
            throw new XMLParseException(e.getMessage());
        }

    }
}

