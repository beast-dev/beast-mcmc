/*
 * PathSamplingAnalysis.java
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

package dr.inference.trace;

import dr.inferencexml.trace.MarginalLikelihoodAnalysisParser;
import dr.util.Attribute;
import dr.util.FileHelpers;
import dr.xml.*;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.util.*;

/**
 * @author Marc A. Suchard
 * @author Alexander Alekseyenko
 * @author Guy Baele
 */
public class PathSamplingAnalysis {

    public static final String PATH_SAMPLING_ANALYSIS = "pathSamplingAnalysis";
    public static final String RESULT_FILE_NAME = "resultsFileName";
    public static final String LIKELIHOOD_COLUMN = "likelihoodColumn";
    public static final String THETA_COLUMN = "thetaColumn";
    public static final String FORMAT = "%5.5g";

    PathSamplingAnalysis(String logLikelihoodName, List<Double> logLikelihoodSample, List<Double> thetaSample) {
        this.logLikelihoodSample = logLikelihoodSample;
        this.logLikelihoodName = logLikelihoodName;
        this.thetaSample = thetaSample;
    }

    public double getLogBayesFactor() {
        if (!logBayesFactorCalculated) {
            calculateBF();
        }
        return logBayesFactor;
    }

    private void calculateBF() {

//  R code from Alex Alekseyenko
//
//  psMLE = function(likelihood, pathParameter){
//      y=tapply(likelihood,pathParameter,mean)
//      L = length(y)
//      midpoints = (y[1:(L-1)] + y[2:L])/2
//      x = as.double(names(y))
//      widths = (x[2:L] - x[1:(L-1)])
//      sum(widths*midpoints)
//  }

        Map<Double, List<Double>> map = new HashMap<Double, List<Double>>();
        orderedTheta = new ArrayList<Double>();

        for (int i = 0; i < logLikelihoodSample.size(); i++) {
            if (!map.containsKey(thetaSample.get(i))) {
                map.put(thetaSample.get(i), new ArrayList<Double>());
                orderedTheta.add(thetaSample.get(i));
            }
            map.get(thetaSample.get(i)).add(logLikelihoodSample.get(i));
        }

        Collections.sort(orderedTheta);

        meanLogLikelihood = new ArrayList<Double>();
        for (double t : orderedTheta) {
            double totalMean = 0;
            int lengthMean = 0;
            List<Double> values = map.get(t);
            for (double v : values) {
                totalMean += v;
                lengthMean++;
            }
            meanLogLikelihood.add(totalMean / lengthMean);
        }

        mlContribution = new ArrayList<Double>();
        
        logBayesFactor = 0;
        innerArea = 0;
        for (int i = 0; i < meanLogLikelihood.size() - 1; i++) {
        	double contribution = (meanLogLikelihood.get(i + 1) + meanLogLikelihood.get(i)) / 2.0 *
            (orderedTheta.get(i + 1) - orderedTheta.get(i));
            logBayesFactor += contribution;
            mlContribution.add(contribution);
            if (i > 0 && i < (meanLogLikelihood.size() - 1)) {
                innerArea += (meanLogLikelihood.get(i + 1) + meanLogLikelihood.get(i)) / 2.0 *
                        (orderedTheta.get(i + 1) - orderedTheta.get(i));
            }
        }
        logBayesFactorCalculated = true;
    }

    public String toString() {
        double bf = getLogBayesFactor();
        StringBuffer sb = new StringBuffer();
        sb.append("PathParameter\tMeanPathLikelihood\tMLContribution\n");
        for (int i = 0; i < orderedTheta.size(); ++i) {
            sb.append(String.format(FORMAT, orderedTheta.get(i)));
            sb.append("\t");
            sb.append(String.format(FORMAT, meanLogLikelihood.get(i)));
            sb.append("\t");
            if (i != (orderedTheta.size()-1)) {
            	sb.append(String.format(FORMAT, mlContribution.get(i)));
            }
            sb.append("\n");
        }

        sb.append("\nlog marginal likelihood (using path sampling) from " + logLikelihoodName + " = " + bf + "\n");
        sb.append("\nInner area for path parameter in ("
                + String.format(FORMAT, orderedTheta.get(1)) + ","
                + String.format(FORMAT, orderedTheta.get(orderedTheta.size() - 2)) + ") = "
                + String.format(FORMAT, innerArea));
        sb.append("\n");
        return sb.toString();
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return PATH_SAMPLING_ANALYSIS;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            String fileName = xo.getStringAttribute(FileHelpers.FILE_NAME);
            String resultFileName = null;
            if (xo.hasAttribute(RESULT_FILE_NAME)) {
                resultFileName = xo.getStringAttribute(RESULT_FILE_NAME);
            }
            StringTokenizer tokenFileName = new StringTokenizer(fileName);
    		int numberOfFiles = tokenFileName.countTokens();
    		System.out.println(numberOfFiles + " file(s) found with marginal likelihood samples");
            try {
            	
            	String likelihoodName = "";
    			List sampleLogLikelihood = null;
    			List sampleTheta = null;
    			
    			for (int j = 0; j < numberOfFiles; j++) {
    				
    				File file = new File(tokenFileName.nextToken());
                    String name = file.getName();
                    String parent = file.getParent();

                    if (!file.isAbsolute()) {
                        parent = System.getProperty("user.dir");
                    }

                    final String fileNamePrefix = System.getProperty("file.name.prefix");
                    final String fileSeparator = System.getProperty("file.separator");
                    if (fileNamePrefix != null) {
                        if (fileNamePrefix.trim().length() == 0 || fileNamePrefix.contains(fileSeparator)) {
                            throw new XMLParseException("The specified file name prefix is illegal.");
                        }
                        file = new File(parent, fileNamePrefix+name);
                    } else {
                        file = new File(parent, name);
                    }

                    fileName = file.getAbsolutePath();
                    
                    XMLObject cxo = xo.getChild(LIKELIHOOD_COLUMN);
                    likelihoodName = cxo.getStringAttribute(Attribute.NAME);
    				
                    cxo = xo.getChild(THETA_COLUMN);
                    String thetaName = cxo.getStringAttribute(Attribute.NAME);
                    
                    LogFileTraces traces = new LogFileTraces(fileName, file);
                    traces.loadTraces();
                    long maxState = traces.getMaxState();
                    
                    // leaving the burnin attribute off will result in 10% being used
                    long burnin = xo.getAttribute(MarginalLikelihoodAnalysisParser.BURN_IN, maxState / 5);

                    if (burnin < 0 || burnin >= maxState) {
                        burnin = maxState / 5;
                        System.out.println("WARNING: Burn-in larger than total number of states - using 20%");
                    }

                    burnin = 0;

                    traces.setBurnIn(burnin);
                    
                    int traceIndexLikelihood = -1;
                    int traceIndexTheta = -1;
                    for (int i = 0; i < traces.getTraceCount(); i++) {
                        String traceName = traces.getTraceName(i);
                        if (traceName.trim().equals(likelihoodName)) {
                            traceIndexLikelihood = i;
                        }
                        if (traceName.trim().equals(thetaName)) {
                            traceIndexTheta = i;
                        }
                    }

                    if (traceIndexLikelihood == -1) {
                        throw new XMLParseException("Column '" + likelihoodName + "' can not be found for " + getParserName() + " element.");
                    }

                    if (traceIndexTheta == -1) {
                        throw new XMLParseException("Column '" + thetaName + "' can not be found for " + getParserName() + " element.");
                    }
                    
                    if (sampleLogLikelihood == null && sampleTheta == null) {
    					sampleLogLikelihood = traces.getValues(traceIndexLikelihood);
    					sampleTheta = traces.getValues(traceIndexTheta);
    				} else {
    					sampleLogLikelihood.addAll(traces.getValues(traceIndexLikelihood));
    					sampleTheta.addAll(traces.getValues(traceIndexTheta));
    				}
                    
    			}

                PathSamplingAnalysis analysis = new PathSamplingAnalysis(likelihoodName, sampleLogLikelihood, sampleTheta);

                System.out.println(analysis.toString());

                if (resultFileName != null) {
                    FileWriter fw = new FileWriter(resultFileName, true);
                    BufferedWriter bw = new BufferedWriter(fw);
                    bw.write(analysis.toString());
                    bw.flush();
                    bw.close();
                }

                return analysis;

            } catch (FileNotFoundException fnfe) {
                throw new XMLParseException("File '" + fileName + "' can not be opened for " + getParserName() + " element.");
            } catch (java.io.IOException ioe) {
                throw new XMLParseException(ioe.getMessage());
            } catch (TraceException e) {
                throw new XMLParseException(e.getMessage());
            }
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "Performs a trace analysis.";
        }

        public Class getReturnType() {
            return PathSamplingAnalysis.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                new StringAttributeRule(FileHelpers.FILE_NAME,
                        "The traceName of a BEAST log file (can not include trees, which should be logged separately)"),
                new StringAttributeRule(RESULT_FILE_NAME,
                        "The name of the output file to which the path sampling estimate will be written", true),
                new ElementRule(THETA_COLUMN, new XMLSyntaxRule[]{
                        new StringAttributeRule(Attribute.NAME, "The column name")}),
                new ElementRule(LIKELIHOOD_COLUMN, new XMLSyntaxRule[]{
                        new StringAttributeRule(Attribute.NAME, "The column name")}),
        };
    };

    private boolean logBayesFactorCalculated = false;
    private double logBayesFactor;
    private double innerArea;
    private final List<Double> logLikelihoodSample;
    private final List<Double> thetaSample;
    private List<Double> meanLogLikelihood;
    private List<Double> mlContribution;
    private final String logLikelihoodName;
    List<Double> orderedTheta;
}
