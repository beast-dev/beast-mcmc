/*
 * TraceAnalysisParser.java
 *
 * Copyright (C) 2002-2006 Alexei Drummond and Andrew Rambaut
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
import dr.util.NumberFormatter;
import dr.xml.*;

import java.io.Reader;
import java.io.File;
import java.io.FileReader;
import java.io.FileNotFoundException;

/**
 *
 * @author Alexei Drummond
 * @version $Id: TraceAnalysisParser.java,v 1.18 2005/05/24 20:26:00 rambaut Exp $
 */
public class TraceAnalysisParser extends AbstractXMLObjectParser {
	
	public static final String TRACE_ANALYSIS = "traceAnalysis";
    public static final String FILE_NAME = "fileName";
	public static final String BURN_IN = "burnIn";
	
	public String getParserName() { return TRACE_ANALYSIS; }
			
	public Object parseXMLObject(XMLObject xo) throws XMLParseException {
			
		try {			
            Reader reader = null;

            String fileName = xo.getStringAttribute(FILE_NAME);
            try {
                File file = new File(fileName);
                String name = file.getName();
                String parent = file.getParent();

                if (!file.isAbsolute()) {
                    parent = System.getProperty("user.dir");
                }

//					System.out.println("Writing log file to "+parent+System.getProperty("path.separator")+name);
                reader = new FileReader(new File(parent, name));
            } catch (FileNotFoundException fnfe) {
                throw new XMLParseException("File '" + fileName + "' can not be opened for " + getParserName() + " element.");
            }

			int burnin = -1;
			if (xo.hasAttribute(BURN_IN)) {
				// leaving the burnin attribute off will result in 10% being used
				burnin = xo.getIntegerAttribute(BURN_IN);
			}
					
			TraceAnalysis[] analysis = TraceAnalysis.report(reader, burnin);
							
			for (int x = 0; x < xo.getChildCount(); x++) {
				XMLObject child = (XMLObject)xo.getChild(x);
				String statName = child.getStringAttribute(Attribute.NAME);
				double expectation = child.getDoubleAttribute(Attribute.VALUE);
				NumberFormatter formatter = new NumberFormatter(6);
				formatter.setPadding(true);
				formatter.setFieldWidth(14);
			
				for (int i =0; i < analysis.length; i++) {
					if (analysis[i].getName().equals(statName)) {
						double estimate = analysis[i].getMean();
						double error = analysis[i].getStdError();
								
						System.out.println("E[" + statName + "]=" + formatter.format(expectation));
									
						if (expectation > (estimate-(2*error)) && expectation < (estimate+(2*error))) {
							System.out.println("OK:       " + formatter.format(estimate) + " +- " + formatter.format(error) + "\n");
						} else {
							System.out.print("WARNING: " + formatter.format(estimate) + " +- " + formatter.format(error) + "\n");
						}
									
					}
				}
			}
					
			System.out.println();
			System.out.flush();
			return analysis;
		} catch (java.io.IOException ioe) {
			throw new XMLParseException(ioe.getMessage());
		}
	}
	
	//************************************************************************
	// AbstractXMLObjectParser implementation
	//************************************************************************
	
	public String getParserDescription() {
		return "Performs a trace analysis. Estimates the mean of the various statistics in the given log file.";
	}
		
	public Class getReturnType() { return TraceAnalysis[].class; }
	
	public XMLSyntaxRule[] getSyntaxRules() { return rules; }
		
	private XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
		new StringAttributeRule(FILE_NAME, "The name of a BEAST log file (can not include trees, which should be logged separately" ),
		AttributeRule.newIntegerRule("burnIn", true)
			//, "The number of states (not sampled states, but actual states) that are discarded from the beginning of the trace before doing the analysis" ),
	};
}