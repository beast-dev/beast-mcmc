/*
 * ArithmeticMeanAnalysisParser.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
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
 *
 */

package dr.inferencexml.trace;

import dr.inference.trace.LogFileTraces;
import dr.inference.trace.MarginalLikelihoodAnalysis;
import dr.inference.trace.TraceException;
import dr.util.Attribute;
import dr.xml.*;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.List;
import java.util.StringTokenizer;

/**
* 
* @author Guy Baele
* 
*/
public class ArithmeticMeanAnalysisParser extends AbstractXMLObjectParser {

	public static final String ARITHMETIC_MEAN_ANALYSIS = "arithmeticMeanAnalysis";
    public static final String FILE_NAME = "fileName";
    public static final String BURN_IN = "burnIn";
    public static final String COLUMN_NAME = "likelihoodColumn";
    public static final String BOOTSTRAP_LENGTH = "bootstrapLength";
    
    public String getParserName() {
        return ARITHMETIC_MEAN_ANALYSIS;
    }
    
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        // Set bootstrapLength
        int bootstrapLength = xo.getAttribute(BOOTSTRAP_LENGTH, 1000);

        long burnin = xo.getAttribute(BURN_IN, 0);

    	String fileName = xo.getStringAttribute(FILE_NAME);
        // Find likelihood column
        XMLObject cxo = xo.getChild(COLUMN_NAME);
        String likelihoodName = cxo.getStringAttribute(Attribute.NAME);

        StringTokenizer tokenFileName = new StringTokenizer(fileName);
        int numberOfFiles = tokenFileName.countTokens();

        ArrayList<Double> sampleLogLikelihood = new ArrayList<Double>();

        // Set analysisType
        String analysisType = "arithmetic";

        try {

            for (int j = 0; j < numberOfFiles; j++) {

                // Open file
                File file = new File(tokenFileName.nextToken());
                String name = file.getName();

                System.out.println("Parsing samples from file: " + name);

                String parent = file.getParent();
                if (!file.isAbsolute()) {
                    parent = System.getProperty("user.dir");
                }
                file = new File(parent, name);
                fileName = file.getAbsolutePath();

                // Load traces and remove burnin
                LogFileTraces traces = new LogFileTraces(fileName, file);
                traces.loadTraces();

                long maxState = traces.getMaxState();

                if (burnin < 0 || burnin >= maxState) {
                    burnin = maxState / 10;
                    System.out.println("WARNING: Burn-in larger than total number of states - using to 10%");
                }

                traces.setBurnIn(burnin);

                int traceIndex = -1;

                for (int i = 0; i < traces.getTraceCount(); i++) {
                    String traceName = traces.getTraceName(i);
                    if (traceName.equals(likelihoodName)) {
                        traceIndex = i;
                        break;
                    }
                }

                if (traceIndex == -1) {
                    throw new XMLParseException("Column '" + likelihoodName + "' can not be found for " + getParserName() + " element.");
                }

                // Get samples and perform analysis
                List<Double> sample = (List)traces.getValues(traceIndex);

                sampleLogLikelihood.addAll(sample);

            }

            System.out.println("Total number of collected samples: " + sampleLogLikelihood.size());

            MarginalLikelihoodAnalysis analysis = new MarginalLikelihoodAnalysis(sampleLogLikelihood,
                    likelihoodName, (int)burnin, analysisType, bootstrapLength);

            System.out.println(analysis.toString());

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
        return "Performs a trace analysis. Estimates the mean of the various statistics in the given log file.";
    }

    public Class getReturnType() {
        return MarginalLikelihoodAnalysis.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new StringAttributeRule(FILE_NAME, "The traceName of a BEAST log file (can not include trees, which should be logged separately"),
            AttributeRule.newIntegerRule(BURN_IN, true),
            AttributeRule.newIntegerRule(BOOTSTRAP_LENGTH, true),
            //, "The number of states (not sampled states, but actual states) that are discarded from the beginning of the trace before doing the analysis" ),
            new ElementRule(COLUMN_NAME, new XMLSyntaxRule[] {
                    new StringAttributeRule(Attribute.NAME,"The column name")}),
    };
	
}
