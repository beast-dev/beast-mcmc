/*
 * GeneralizedSteppingStoneSamplingAnalysis.java
 *
 * Copyright (c) 2002-2021 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

import dr.util.Attribute;
import dr.util.FileHelpers;
import dr.xml.*;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.nio.file.Paths;
import java.util.*;

/**
 * @author Guy Baele
 */

public class GeneralizedHarmonicMeanAnalysis implements Reportable {

    public static final String GENERALIZED_HARMONIC_MEAN_ANALYSIS = "generalizedHarmonicMeanAnalysis";
    public static final String RESULT_FILE_NAME = "resultsFileName";
    public static final String SOURCE_COLUMN = "sourceColumn";
    public static final String DESTINATION_COLUMN = "destinationColumn";
    public static final String FORMAT = "%5.5g";

    private final String sourceName, destinationName;
    private final List<Double> sourceSample;
    private final List<Double> destinationSample;
    private boolean logBayesFactorCalculated = false;
    private double logBayesFactor;
    private List<Double> maxLogLikelihood;
    private List<Double> mlContribution;

    public GeneralizedHarmonicMeanAnalysis(String sourceName, String destinationName, List<Double> sourceSample, List<Double> destinationSample) {
        this.sourceName = sourceName;
        this.destinationName = destinationName;
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
        
        logBayesFactor = 0.0;

        if (this.sourceSample.size() != this.destinationSample.size()) {
            throw new RuntimeException("Error in GeneralizedHarmonicMeanAnalysis: source and destination samples have different numbers of entries.");
        }

        int size = this.sourceSample.size();

        ArrayList<Double> exponentialParts = new ArrayList<Double>();
        for (int i = 0; i < size; i++) {
            double exponentialPart = destinationSample.get(i)*(1.0 - sourceSample.get(i)/destinationSample.get(i));
            exponentialParts.add(exponentialPart);
            //System.out.println("exponentialPart = " + exponentialPart);
        }
        double firstPart = exponentialParts.get(0);

        double internalSum = 1.0;
        for (int i = 1; i < exponentialParts.size(); i++) {
            internalSum += Math.exp(exponentialParts.get(i) - firstPart);
            System.out.println("internalSum = " + internalSum);
        }

        System.out.println("\nlog(size) = " + StrictMath.log(size));
        System.out.println("-firstPart = " +  (-firstPart));
        System.out.println("-Math.log(internalSum) = " + (-Math.log(internalSum)));

        logBayesFactor = StrictMath.log(size) - firstPart - Math.log(internalSum);

        logBayesFactorCalculated = true;
            
    }
    
    public String toString() {
        double bf = getLogBayesFactor();
        StringBuffer sb = new StringBuffer();
        sb.append("\nlog marginal likelihood (using generalized harmonic mean) from (" + sourceName + " - " + destinationName + ") = " + bf + "\n");
        return sb.toString();
    }

    @Override
    public String getReport() {
        return this.toString();
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {
    	
    	public String getParserName() {
            return GENERALIZED_HARMONIC_MEAN_ANALYSIS;
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
    			List sampleSource = null;
    			List sampleDestination = null;
    			
    			for (int j = 0; j < numberOfFiles; j++) {
    			
    				File file = new File(tokenFileName.nextToken());
    				String name = file.getName();
    				String parent = file.getParent();

                    if (!file.isAbsolute()) {
                        if (parent == null) {
                            parent = System.getProperty("user.dir");
                        } else {
                            parent = Paths.get(System.getProperty("user.dir"), parent).toString();
                        }
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

    				LogFileTraces traces = new LogFileTraces(fileName, file);
    				traces.loadTraces();
                
    				long burnin = 0;
                
    				traces.setBurnIn(burnin);

    				int traceIndexSource = -1;
    				int traceIndexDestination = -1;
    				for (int i = 0; i < traces.getTraceCount(); i++) {
    					String traceName = traces.getTraceName(i);
    					if (traceName.trim().equals(sourceName)) {
    						traceIndexSource = i;
    					}
    					if (traceName.trim().equals(destinationName)) {
    						traceIndexDestination = i;
    					}
    				}

    				if (traceIndexSource == -1) {
    					throw new XMLParseException("Column '" + sourceName + "' can not be found for " + getParserName() + " element.");
    				}
    				
    				if (traceIndexDestination == -1) {
    					throw new XMLParseException("Column '" + destinationName + "' can not be found for " + getParserName() + " element.");
    				}

    				if (sampleSource == null && sampleDestination == null) {
    					sampleSource = traces.getValues(traceIndexSource);
    					sampleDestination = traces.getValues(traceIndexDestination);
    				} else {
    					sampleSource.addAll(traces.getValues(traceIndexSource));
    					sampleDestination.addAll(traces.getValues(traceIndexDestination));
    				}
    			
    			}
    			
                GeneralizedHarmonicMeanAnalysis analysis = new GeneralizedHarmonicMeanAnalysis(sourceName, destinationName, sampleSource, sampleDestination);

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
            return GeneralizedHarmonicMeanAnalysis.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                new StringAttributeRule(FileHelpers.FILE_NAME,
                        "The traceName of a BEAST log file (can not include trees, which should be logged separately)"),
                new StringAttributeRule(RESULT_FILE_NAME,
                        "The name of the output file to which the generalized stepping-stone sampling estimate will be written", true),
                /*new ElementRule(THETA_COLUMN, new XMLSyntaxRule[]{
                        new StringAttributeRule(Attribute.NAME, "The column name")}),*/
                new ElementRule(SOURCE_COLUMN, new XMLSyntaxRule[]{
                        new StringAttributeRule(Attribute.NAME, "The column name")}),
                new ElementRule(DESTINATION_COLUMN, new XMLSyntaxRule[]{
                        new StringAttributeRule(Attribute.NAME, "The column name")})
        };
    	
    };
	
}