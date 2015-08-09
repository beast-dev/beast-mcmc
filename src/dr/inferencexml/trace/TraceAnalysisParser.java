/*
 * TraceAnalysisParser.java
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

import dr.inference.trace.*;
import dr.util.Attribute;
import dr.util.NumberFormatter;
import dr.xml.*;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.List;

/**
 * @author Alexei Drummond
 * @version $Id: TraceAnalysisParser.java,v 1.18 2005/05/24 20:26:00 rambaut Exp $
 */
public class TraceAnalysisParser extends AbstractXMLObjectParser {

    public static final String TRACE_ANALYSIS = "traceAnalysis";
    public static final String FILE_NAME = "fileName";
    public static final String BURN_IN = "burnIn";
    public static final String STD_ERROR = "stdError";
    public static final String EXPECTATION = "expectation";
    public static final String COMPUTE_MSE = "computeMSE";

    public String getParserName() {
        return TRACE_ANALYSIS;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        String fileName = xo.getStringAttribute(FILE_NAME);
        boolean withStdError = xo.getAttribute(STD_ERROR, false);
        boolean computeMSE = xo.getAttribute(COMPUTE_MSE, false);
        try {

            File file = new File(fileName);
            String name = file.getName();
            String parent = file.getParent();

            if (!file.isAbsolute()) {
                parent = System.getProperty("user.dir");
            }

            file = new File(parent + File.separator, name);
            if (file.exists()) {
                fileName = file.getName();

                // leaving the burnin attribute off will result in 10% being used
                int burnin = xo.getAttribute(BURN_IN, -1);

                TraceList traces = TraceAnalysis.report(fileName, burnin, null, withStdError);
                for (int x = 0; x < xo.getChildCount(); x++) {
                    XMLObject child = (XMLObject) xo.getChild(x);
                    String statName = child.getStringAttribute(Attribute.NAME);
                    double expectation = child.getDoubleAttribute(Attribute.VALUE);
                    NumberFormatter formatter = new NumberFormatter(6);
                    formatter.setPadding(true);
                    formatter.setFieldWidth(14);

                    for (int i = 0; i < traces.getTraceCount(); i++) {
                        TraceDistribution distribution = traces.getDistributionStatistics(i);
                        TraceCorrelation corr = traces.getCorrelationStatistics(i);
                        if (traces.getTraceName(i).equals(statName)) {
                            double estimate = distribution.getMean();
                            double error = corr.getStdErrorOfMean();

                            System.out.print("E[" + statName + "] = " + formatter.format(expectation));

                            if (computeMSE) {
                                List values = traces.getValues(i);
                                double[] dv = new double[values.size()];
                                for (int j = 0; j < dv.length; j++) {
                                    dv[j] = ((Number) values.get(j)).doubleValue();
                                }
                                double MSE = distribution.getMeanSquaredError(dv, expectation);
                                System.out.println(" MSE = " + formatter.format(MSE));
                            } else {
                                System.out.println("");


                                if (expectation > (estimate - (2 * error)) && expectation < (estimate + (2 * error))) {
                                    System.out.println("OK:       " + formatter.format(estimate) + " +- " + formatter.format(error) + "\n");
                                } else {
                                    System.out.print("WARNING: " + formatter.format(estimate) + " +- " + formatter.format(error) + "\n");
                                }
                            }
                        }
                    }
                }

                System.out.println();
                System.out.flush();
                return traces;
            } else {
                throw new XMLParseException("Log file, " + parent + File.separator + name + " does not exist.");
            }
        } catch (FileNotFoundException fnfe) {
            throw new XMLParseException("File '" + fileName + "' can not be opened for " + getParserName() + " element.");
        } catch (java.io.IOException ioe) {
            throw new XMLParseException(ioe.getMessage());
        } catch (TraceException e) {
            throw new XMLParseException(e.toString());
        }
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "Performs a trace analysis. Estimates the mean of the various statistics in the given log file.";
    }

    public Class getReturnType() {
        return TraceAnalysis[].class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new StringAttributeRule(FILE_NAME, "The name of a BEAST log file (can not include trees, which should be logged separately"),
            AttributeRule.newIntegerRule(BURN_IN, true),
            new ElementRule(EXPECTATION, new XMLSyntaxRule[]{AttributeRule.newStringRule("name"), AttributeRule.newStringRule("value")}, 0, Integer.MAX_VALUE),
            AttributeRule.newBooleanRule(STD_ERROR, true),
            AttributeRule.newBooleanRule(COMPUTE_MSE, true),
    };
}