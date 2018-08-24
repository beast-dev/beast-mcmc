/*
 * ArithmeticMeanAnalysisParser.java
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

package dr.inferencexml.trace;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.List;

import dr.inference.trace.LogFileTraces;
import dr.inference.trace.MarginalLikelihoodAnalysis;
import dr.inference.trace.TraceException;
import dr.util.Attribute;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.AttributeRule;
import dr.xml.ElementRule;
import dr.xml.StringAttributeRule;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

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
    	
    	String fileName = xo.getStringAttribute(FILE_NAME);
        try {

            // Open file
            File file = new File(fileName);
            String name = file.getName();
            String parent = file.getParent();
            if (!file.isAbsolute()) {
                parent = System.getProperty("user.dir");
            }
            file = new File(parent, name);
            fileName = file.getAbsolutePath();

            // Set analysisType
            String analysisType = "arithmetic";
            
            // Set bootstrapLength
            int bootstrapLength = xo.getAttribute(BOOTSTRAP_LENGTH, 1000);
            
            // Load traces and remove burnin
            LogFileTraces traces = new LogFileTraces(fileName, file);
            traces.loadTraces();
            long maxState = traces.getMaxState();

            long burnin = xo.getAttribute(BURN_IN, maxState / 10);

            if (burnin < 0 || burnin >= maxState) {
                burnin = maxState / 10;
                System.out.println("WARNING: Burn-in larger than total number of states - using to 10%");
            }

            traces.setBurnIn(burnin);
            
            // Find likelihood column
            XMLObject cxo = xo.getChild(COLUMN_NAME);
            String likelihoodName = cxo.getStringAttribute(Attribute.NAME);
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

            MarginalLikelihoodAnalysis analysis = new MarginalLikelihoodAnalysis(sample,
                    traces.getTraceName(traceIndex), (int)burnin, analysisType, bootstrapLength);

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
