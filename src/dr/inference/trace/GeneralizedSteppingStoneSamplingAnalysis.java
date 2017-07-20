/*
 * GeneralizedSteppingStoneSamplingAnalysis.java
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

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.util.*;

import dr.util.Attribute;
import dr.util.FileHelpers;

import dr.xml.*;

/**
 * @author Guy Baele
 */

public class GeneralizedSteppingStoneSamplingAnalysis {
	
    public static final String GENERALIZED_STEPPING_STONE_SAMPLING_ANALYSIS = "generalizedSteppingStoneSamplingAnalysis";
    public static final String RESULT_FILE_NAME = "resultsFileName";
    public static final String THETA_COLUMN = "thetaColumn";
    public static final String SOURCE_COLUMN = "sourceColumn";
    public static final String DESTINATION_COLUMN = "destinationColumn";
    public static final String FORMAT = "%5.5g";
    
    private final String sourceName, destinationName;
    private final List<Double> thetaSample;
    private final List<Double> sourceSample;
    private final List<Double> destinationSample;
    private boolean logBayesFactorCalculated = false;
    private double logBayesFactor;
    private List<Double> maxLogLikelihood;
    private List<Double> orderedTheta;
    private List<Double> mlContribution;
    
    public GeneralizedSteppingStoneSamplingAnalysis(String sourceName, String destinationName, List<Double> thetaSample, List<Double> sourceSample, List<Double> destinationSample) {
        this.sourceName = sourceName;
        this.destinationName = destinationName;
        this.thetaSample = thetaSample;
        this.sourceSample = sourceSample;
        this.destinationSample = destinationSample;
    }
    
    public double getLogBayesFactor() {
        if (!logBayesFactorCalculated) {
            calculateBF();
        }
        return logBayesFactor;
    }
    
    private void calculateBF() {
    	
    	Map<Double, List<Double>> map = new HashMap<Double, List<Double>>();
    	Map<Double, List<Double>> testmap = new HashMap<Double, List<Double>>();
        orderedTheta = new ArrayList<Double>();
        
        //the log-likelihood*prior/refprior values are needed to calculate the marginal likelihood
        for (int i = 0; i < sourceSample.size(); i++) {
        	if (!map.containsKey(thetaSample.get(i))) {
        		map.put(thetaSample.get(i), new ArrayList<Double>());
        		testmap.put(thetaSample.get(i), new ArrayList<Double>());
                orderedTheta.add(thetaSample.get(i));
            }
            map.get(thetaSample.get(i)).add(sourceSample.get(i) - destinationSample.get(i));
            testmap.get(thetaSample.get(i)).add(sourceSample.get(i));
        }

        Collections.sort(orderedTheta);

        //a list with the maxima of the log-likelihood*prior/refprior values is constructed
        System.out.println("Test source column:");
        maxLogLikelihood = new ArrayList<Double>();
        for (double t : orderedTheta) {
        	List<Double> values = map.get(t);
            maxLogLikelihood.add(Collections.max(values));  
            System.out.println(Collections.max(testmap.get(t)));
        }
        
        System.out.println("Number of maximum loglikelihoods: " + maxLogLikelihood.size());
        for (double ml : maxLogLikelihood) {
        	System.out.println(ml);
        }

        mlContribution = new ArrayList<Double>();
        
        logBayesFactor = 0.0;
        for (int i = 1; i < orderedTheta.size(); i++) {
        	double contribution = (orderedTheta.get(i) - orderedTheta.get(i-1)) * maxLogLikelihood.get(i-1);
        	logBayesFactor += contribution;
        	mlContribution.add(contribution);
        }
        //System.out.println(logBayesFactor);
        
        for (int i = 1; i < orderedTheta.size(); i++) {
        	double internalSum = 0.0;
        	for (int j = 0; j < map.get(orderedTheta.get(i-1)).size(); j++) {
        		internalSum += Math.exp((orderedTheta.get(i) - orderedTheta.get(i-1)) * (map.get(orderedTheta.get(i-1)).get(j) - maxLogLikelihood.get(i-1)));
        	}
        	internalSum /= map.get(orderedTheta.get(i-1)).size();
        	//System.out.print(orderedTheta.get(i) + "-" + orderedTheta.get(i-1) + ": ");
        	//System.out.println(Math.log(internalSum));
        	mlContribution.set(i-1, mlContribution.get(i-1) + Math.log(internalSum));
        	logBayesFactor += Math.log(internalSum);
        }
        
        logBayesFactorCalculated = true;
            
    }
    
    public String toString() {
        double bf = getLogBayesFactor();
        StringBuffer sb = new StringBuffer();
        sb.append("PathParameter\tMaxPathLikelihood\tMLContribution\n");
        for (int i = 0; i < orderedTheta.size(); ++i) {
            sb.append(String.format(FORMAT, orderedTheta.get(i)));
            sb.append("\t");
            sb.append(String.format(FORMAT, maxLogLikelihood.get(i)));
            sb.append("\t");
            if (i != (orderedTheta.size()-1)) {
            	sb.append(String.format(FORMAT, mlContribution.get(i)));
            }
            sb.append("\n");
        }

        sb.append("\nlog marginal likelihood (using generalized stepping stone sampling) from (" + sourceName + " - " + destinationName + ") = " + bf + "\n");
        return sb.toString();
    }
    
    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {
    	
    	public String getParserName() {
            return GENERALIZED_STEPPING_STONE_SAMPLING_ANALYSIS;
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
    			
    			String sourceName = "", destinationName = "";
    			List sampleTheta = null;
    			List sampleSource = null;
    			List sampleDestination = null;
    			
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

    				XMLObject cxo = xo.getChild(SOURCE_COLUMN);
    				sourceName = cxo.getStringAttribute(Attribute.NAME);
    				
    				cxo = xo.getChild(DESTINATION_COLUMN);
    				destinationName = cxo.getStringAttribute(Attribute.NAME);

    				cxo = xo.getChild(THETA_COLUMN);
    				String thetaName = cxo.getStringAttribute(Attribute.NAME);

    				LogFileTraces traces = new LogFileTraces(fileName, file);
    				traces.loadTraces();
                
    				long burnin = 0;
                
    				traces.setBurnIn(burnin);

    				int traceIndexTheta = -1;
    				int traceIndexSource = -1;
    				int traceIndexDestination = -1;
    				for (int i = 0; i < traces.getTraceCount(); i++) {
    					String traceName = traces.getTraceName(i);
    					if (traceName.trim().equals(thetaName)) {
    						traceIndexTheta = i;
    					}
    					if (traceName.trim().equals(sourceName)) {
    						traceIndexSource = i;
    					}
    					if (traceName.trim().equals(destinationName)) {
    						traceIndexDestination = i;
    					}
    				}

    				if (traceIndexTheta == -1) {
    					throw new XMLParseException("Column '" + thetaName + "' can not be found for " + getParserName() + " element.");
    				}
    				
    				if (traceIndexSource == -1) {
    					throw new XMLParseException("Column '" + sourceName + "' can not be found for " + getParserName() + " element.");
    				}
    				
    				if (traceIndexDestination == -1) {
    					throw new XMLParseException("Column '" + destinationName + "' can not be found for " + getParserName() + " element.");
    				}

    				if (sampleTheta == null && sampleSource == null && sampleDestination == null) {
    					sampleTheta = traces.getValues(traceIndexTheta);
    					sampleSource = traces.getValues(traceIndexSource);
    					sampleDestination = traces.getValues(traceIndexDestination);
    				} else {
    					sampleTheta.addAll(traces.getValues(traceIndexTheta));
    					sampleSource.addAll(traces.getValues(traceIndexSource));
    					sampleDestination.addAll(traces.getValues(traceIndexDestination));
    				}
    			
    			}
    			
                GeneralizedSteppingStoneSamplingAnalysis analysis = new GeneralizedSteppingStoneSamplingAnalysis(sourceName, destinationName, sampleTheta, sampleSource, sampleDestination);

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
            return GeneralizedSteppingStoneSamplingAnalysis.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                new StringAttributeRule(FileHelpers.FILE_NAME,
                        "The traceName of a BEAST log file (can not include trees, which should be logged separately)"),
                new StringAttributeRule(RESULT_FILE_NAME,
                        "The name of the output file to which the generalized stepping-stone sampling estimate will be written", true),
                new ElementRule(THETA_COLUMN, new XMLSyntaxRule[]{
                        new StringAttributeRule(Attribute.NAME, "The column name")}),
                new ElementRule(SOURCE_COLUMN, new XMLSyntaxRule[]{
                        new StringAttributeRule(Attribute.NAME, "The column name")}),
                new ElementRule(DESTINATION_COLUMN, new XMLSyntaxRule[]{
                        new StringAttributeRule(Attribute.NAME, "The column name")})
        };
    	
    };
	
}