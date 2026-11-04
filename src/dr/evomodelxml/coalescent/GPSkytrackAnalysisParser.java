/*
 * GPSkytrackAnalysisParser.java
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

package dr.evomodelxml.coalescent;

import dr.evolution.io.Importer;
import dr.evomodel.coalescent.GPSkytrackAnalysis;
import dr.inference.model.Parameter;
//import dr.evomodel.coalescent.GaussianProcessSkytrackLikelihood;
import dr.inference.trace.TraceException;
import dr.util.FileHelpers;
import dr.xml.*;

import java.io.File;
//import java.io.FileWriter;
//import java.io.PrintWriter;

/**
 */
public class GPSkytrackAnalysisParser extends AbstractXMLObjectParser {

    public static final String GP_ANALYSIS = "GPAnalysis";
    public static final String FILE_NAME = "fileName";
    public static final String BURN_IN = "burnIn";
//    public static final String HPD_LEVELS = "Confidencelevels";
//    public static final String QUANTILES = "useQuantiles";
//    public static final String LOG_SPACE = VariableDemographicModelParser.LOG_SPACE;
    public static final String N_GRID = "numGridPoints";
//    public static final String N_CHANGES = "nChanges";

//    public static final String TREE_LOG = "treeOfLoci";

    public static final String LOG_FILE_NAME = "logFileName";

    private String getElementText(XMLObject xo, String childName) throws XMLParseException {
        return xo.getChild(childName).getStringChild(0);
    }

    public String getParserName() {
        return GP_ANALYSIS;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        System.err.println("The Summary Statistics are being created...");

        try {

            // 10% is brun-in default
            final double burnin = xo.getAttribute(BURN_IN, 0.1);
            if (burnin < 0) {
                throw new XMLParseException("burnIn should be either between 0 and 1 or a positive number");
            }



            Parameter numGridPoints = new Parameter.Default(0,1);
            if (xo.getChild(N_GRID) != null) {
                XMLObject cxo = xo.getChild(N_GRID);
                numGridPoints = (Parameter) cxo.getChild(Parameter.class);
            }


            final File log = FileHelpers.getFile(getElementText(xo, LOG_FILE_NAME));


            return new dr.evomodel.coalescent.GPSkytrackAnalysis(log,burnin,numGridPoints);

        } catch (java.io.IOException ioe) {
            throw new XMLParseException(ioe.getMessage());
        } catch (Importer.ImportException e) {
            throw new XMLParseException(e.toString());
        } catch (TraceException e) {
            throw new XMLParseException(e.toString());
        }
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "reconstruct population graph from GPSkytrack run.";
    }

    public Class getReturnType() {
        return GPSkytrackAnalysis.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newDoubleRule(BURN_IN, true, "The number of states (not sampled states, but" +
                    " actual states) that are discarded from the beginning of the trace and are excluded from " +
                    "the analysis"),
            new ElementRule(LOG_FILE_NAME, String.class, "The name of a BEAST log file"),
            new ElementRule(N_GRID, new XMLSyntaxRule[]{
             new ElementRule(Parameter.class)
    })
    };

}
