/*
 * CnCsPerSiteAnalysis.java
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

import dr.util.*;
import dr.xml.*;
import mpi.Comm;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.*;

/**
 * @author Philippe Lemey
 * @author Marc A. Suchard
 */
public class CnCsPerSiteAnalysis implements Citable {
    public static final String CNCS_PER_SITE_ANALYSIS = "cNcSPerSiteAnalysis";
    public static final String BURN_IN = "burnin";
    public static final String CUTOFF = "cutoff";
    public static final String INCLUDE_SIGNIFICANT_SYMBOL = "includeSymbol";
    public static final String INCLUDE_PVALUES = "includePValues";
    public static final String INCLUDE_SITE_CLASSIFICATION = "includeClassification";
    public static final String SEPARATOR_STRING = "separator";
    public static final String INCLUDE_SIMULATION_OUTCOME = "simulationOutcome";
    public static final String SITE_SIMULATION = "siteSimulation";

    public CnCsPerSiteAnalysis(TraceList traceListN, TraceList traceListS) {
        this.traceListN = traceListN;
        this.traceListS = traceListS;
        this.numSites = (traceListN.getTraceCount())/2;
        this.format = new OutputFormat();


        fieldWidth = 14;
        firstField = 10;
        numberFormatter = new NumberFormatter(6);
        numberFormatter.setPadding(true);
        numberFormatter.setFieldWidth(fieldWidth);
    }

    public void setIncludeMeans(boolean b) {
        format.includeMeans = b;
    }

    public void setIncludeSignificantSymbol(boolean b) {
        format.includeSignificantSymbol = b;
    }

    public void setincludePValues(boolean b) {
        format.includePValues = b;
    }

    public void setIncludeSimulationOutcome(boolean b) {
        format.includeSimulationOutcome = b;
    }

    public boolean getIncludeSimulationOutcome() {
        return(format.includeSimulationOutcome);
    }

    public void setSiteSimulation(String[] d) {
        format.siteSimulation = d;
    }

    public void setIncludeSiteClassification(boolean b) {
        format.includeSiteClassification = b;
    }

    public void setCutoff(double d) {
        format.cutoff = d;
    }

    public void setSeparator(String s) {
        format.separator = s;
    }

    private String toStringSite(int index, OutputFormat format) {
        StringBuilder sb = new StringBuilder();
        traceListN.analyseTrace(index);
        traceListN.analyseTrace((traceListN.getTraceCount()/2) + index);
        traceListS.analyseTrace(index);
        traceListS.analyseTrace((traceListS.getTraceCount()/2) + index);
        TraceDistribution distributionCN = traceListN.getDistributionStatistics(index);
        TraceDistribution distributionCS = traceListS.getDistributionStatistics(index);
        TraceDistribution distributionUN = traceListN.getDistributionStatistics((traceListN.getTraceCount()/2) + index);
        TraceDistribution distributionUS = traceListS.getDistributionStatistics((traceListS.getTraceCount()/2) + index);
        double meanCN = distributionCN.getMean();
        double meanCS = distributionCS.getMean();
        double meanUN = distributionUN.getMean();
        double meanUS = distributionUS.getMean();

        sb.append(numberFormatter.formatToFieldWidth(Integer.toString(index + 1), firstField));

        if (format.includeMeans) {
            sb.append(format.separator);
            sb.append(meanCN);
            sb.append(format.separator);
            sb.append(meanUN);
            sb.append(format.separator);
            sb.append(meanCS);
            sb.append(format.separator);
            sb.append(meanUS);
        }

        boolean isSignificant = false;
        String classification = "0";

        double negativeProb =  0;
        if((meanCN+meanCS) == 0){
            negativeProb = 1;
        } else {
            negativeProb =  getCumExtBinProb(meanCS,(meanCN+meanCS),(meanUS/(meanUN+meanUS)));
        }

        if (negativeProb < format.cutoff) {
            isSignificant = true;
            classification = "-";
        }

        double positiveProb =  0;
        if((meanCN+meanCS) == 0){
            positiveProb = 1;
        } else {
            positiveProb =  getCumExtBinProb(meanCN,(meanCN+meanCS),(meanUN/(meanUN+meanUS)));
        }

        if (positiveProb < format.cutoff) {
            isSignificant = true;
            classification = "+";
        }

        if (format.includePValues) {
            sb.append(format.separator);
            sb.append(negativeProb);
            sb.append(format.separator);
            sb.append(positiveProb);
        }

        if (format.includeSiteClassification) {
            sb.append(format.separator);
            sb.append(classification);
        }


        if (format.includeSignificantSymbol) {
            sb.append(format.separator);
            if (isSignificant) {
                sb.append("*");
            } else {
                // Do nothing?
            }
        }

        if (format.includeSimulationOutcome) {
            sb.append(format.separator);
            sb.append(format.siteSimulation[index]);
            sb.append(format.separator);
            if (format.siteSimulation[index].equals("+") || format.siteSimulation[index].equals("-")) {
                if (classification.equals(format.siteSimulation[index])){
                    sb.append("TP");   // True Positive
                } else {
                    sb.append("FN");   // True Negative
                }
            }  else {
                if (classification.equals(format.siteSimulation[index])){
                    sb.append("TN");   // True Negative
                } else {
                    sb.append("FP");   // False Positive
                }
            }
        }


        sb.append("\n");
        return sb.toString();
    }

    public String header(OutputFormat format) {
        StringBuilder sb = new StringBuilder();
        sb.append("# Some information here\n");
        sb.append("# Please cite: " + Utils.getCitationString(this));


        sb.append(numberFormatter.formatToFieldWidth("Site", firstField));

        if (format.includeMeans) {
            sb.append(format.separator);
            sb.append(numberFormatter.formatToFieldWidth("Mean CN", fieldWidth));
            sb.append(numberFormatter.formatToFieldWidth("Mean UN", fieldWidth));
            sb.append(numberFormatter.formatToFieldWidth("Mean CS", fieldWidth));
            sb.append(numberFormatter.formatToFieldWidth("Mean US", fieldWidth));
        }

        if (format.includePValues) {
            sb.append(format.separator);
            sb.append(numberFormatter.formatToFieldWidth("Neg_Pvalue", fieldWidth));
            sb.append(format.separator);
            sb.append(numberFormatter.formatToFieldWidth("Pos_Pvalue", fieldWidth));
        }

        if (format.includeSiteClassification) {
            sb.append(format.separator);
            sb.append(numberFormatter.formatToFieldWidth("Classification", fieldWidth));
        }

        if (format.includeSignificantSymbol) {
            sb.append(format.separator);
            sb.append(numberFormatter.formatToFieldWidth("Significant", fieldWidth));
        }
        if (format.includeSimulationOutcome) {
            sb.append(format.separator);
            sb.append(numberFormatter.formatToFieldWidth("Simulated", fieldWidth));
            sb.append(format.separator);
            sb.append(numberFormatter.formatToFieldWidth("Evaluation", fieldWidth));
        }
        sb.append("\n");
        return sb.toString();
    }

    public String toString() {

        StringBuilder sb = new StringBuilder();
        sb.append(header(format));
        for (int i = 0; i < numSites; ++i) {
            sb.append(toStringSite(i, format));
        }

        return sb.toString();
    }

    @Override
    public Citation.Category getCategory() {
        return Citation.Category.COUNTING_PROCESSES;
    }

    @Override
    public String getDescription() {
        return "Renaissance counting";
    }

    @Override
    public List<Citation> getCitations() {
        return Collections.singletonList(CommonCitations.LEMEY_2012_RENAISSANCE);
    }

    private class OutputFormat {
        boolean includeMeans;
        //        boolean includeHPD;
//        boolean includeSignificanceLevel;
        boolean includePValues;
        boolean includeSignificantSymbol;
        boolean includeSiteClassification;
        boolean includeSimulationOutcome;
        String[] siteSimulation;
        double cutoff;
        //        double proportion;
//        SignificanceTest test;
        String separator;

        OutputFormat() {
            this(true, true, true, true, false, null, 0.05, "\t");
        }

        OutputFormat(boolean includeMeans,
                     boolean includePValues,
                     boolean includeSignificantSymbol,
                     boolean includeSiteClassification,
                     boolean includeSimulationOutcome,
                     String[] siteSimulation,
                     double cutoff,
                     String separator) {
            this.includeMeans = includeMeans;
            this.includePValues = includePValues;
            this.includeSignificantSymbol = includeSignificantSymbol;
            this.includeSiteClassification = includeSiteClassification;
            this.includeSimulationOutcome = includeSimulationOutcome;
            this.siteSimulation = siteSimulation;
            this.cutoff = cutoff;
            this.separator = separator;
        }
    }

    private static double getCumExtBinProb(double x, double N, double P){
        double cumProb = 0;
        return cumProb;
    }

    private static double[] toPrimitiveDoubleArray(Double[] array){
        double[] returnArray = new double[array.length];
        for(int i = 0; i < array.length; i++ ){
            returnArray[i] = array[i].doubleValue();
        }
        return returnArray;
    }

    private static String[] parseVariableLengthStringArray(String inString) {

        List<String> returnList = new ArrayList<String>();
        StringTokenizer st = new StringTokenizer(inString, ",");
        while (st.hasMoreTokens()) {
            returnList.add(st.nextToken());
        }

        if (returnList.size() > 0) {
            String[] stringArray = new String[returnList.size()];
            stringArray = returnList.toArray(stringArray);
            return stringArray;
        }
        return null;
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return CNCS_PER_SITE_ANALYSIS;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            String fileNameCN = xo.getStringAttribute(FileHelpers.FILE_NAME +"CN");
            String fileNameCS = xo.getStringAttribute(FileHelpers.FILE_NAME +"CS");
            try {

                File fileCN = new File(fileNameCN);
                File fileCS = new File(fileNameCS);
                String nameCN = fileCN.getName();
                String nameCS = fileCS.getName();
                String parentCN = fileCN.getParent();
                String parentCS = fileCS.getParent();

                if (!fileCN.isAbsolute()) {
                    parentCN = System.getProperty("user.dir");
                }
                if (!fileCS.isAbsolute()) {
                    parentCS = System.getProperty("user.dir");
                }

                fileCN = new File(parentCN, nameCN);
                fileCS = new File(parentCS, nameCS);

                fileNameCS = fileCS.getAbsolutePath();

                LogFileTraces tracesCN = new LogFileTraces(fileNameCN, fileCN);
                LogFileTraces tracesCS = new LogFileTraces(fileNameCS, fileCS);
                tracesCN.loadTraces();
                tracesCS.loadTraces();

                long maxStateCN = tracesCN.getMaxState();
                long maxStateCS = tracesCS.getMaxState();

                if (maxStateCN != maxStateCS){
                    System.err.println("max states in" + fileNameCN + "and" + fileNameCS + "are not equal");
                }

                // leaving the burnin attribute off will result in 10% being used
                long burnin = xo.getAttribute(BURN_IN, maxStateCN / 10);
                //TODO: implement custom burn-in

                if (burnin < 0 || burnin >= maxStateCN) {
                    burnin = maxStateCN / 5;
                    System.out.println("WARNING: Burn-in larger than total number of states - using 20%");
                }

                tracesCN.setBurnIn(burnin);
                tracesCS.setBurnIn(burnin);

                // TODO: Filter traces to include only dNdS columns

                CnCsPerSiteAnalysis analysis = new CnCsPerSiteAnalysis(tracesCN,tracesCS);

                analysis.setCutoff(xo.getAttribute(CUTOFF, 0.05));
//                analysis.setProportion(xo.getAttribute(PROPORTION, 0.95));
                analysis.setSeparator(xo.getAttribute(SEPARATOR_STRING, "\t"));
                analysis.setincludePValues(xo.getAttribute(INCLUDE_PVALUES, true));
                analysis.setIncludeSignificantSymbol(xo.getAttribute(INCLUDE_SIGNIFICANT_SYMBOL, true));
                analysis.setIncludeSiteClassification(xo.getAttribute(INCLUDE_SITE_CLASSIFICATION, true));
                analysis.setIncludeSimulationOutcome(xo.getAttribute(INCLUDE_SIMULATION_OUTCOME, false));
                if (analysis.getIncludeSimulationOutcome()){
                    String sites = (String)xo.getAttribute(SITE_SIMULATION, "empty");
                    if (sites.equals("empty")){
                        System.err.println("you want simulation evaluation but do not provide a site simulation string??");
                    } else {
                        String[] siteSimulation = parseVariableLengthStringArray(sites);
                        analysis.setSiteSimulation(siteSimulation);
                    }
                }


                return analysis;

            } catch (FileNotFoundException fnfe) {
                throw new XMLParseException("File '" + fileNameCN + "or" + fileNameCS + "' can not be opened for " + getParserName() + " element.");
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
            return "Performs a trace analysis of N and S counts.";
        }

        public Class getReturnType() {
            return CnCsPerSiteAnalysis.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                AttributeRule.newDoubleRule(CUTOFF, true),
                AttributeRule.newIntegerRule(BURN_IN,true),
                AttributeRule.newBooleanRule(INCLUDE_PVALUES, true),
                AttributeRule.newBooleanRule(INCLUDE_SIGNIFICANT_SYMBOL, true),
                AttributeRule.newBooleanRule(INCLUDE_SITE_CLASSIFICATION, true),
                AttributeRule.newBooleanRule(INCLUDE_SIMULATION_OUTCOME, true),
                AttributeRule.newStringRule(SITE_SIMULATION, true),
                AttributeRule.newStringRule(SEPARATOR_STRING, true),
//                new StringAttributeRule(FileHelpers.FILE_NAME,
//                        "The traceName of a BEAST log file (can not include trees, which should be logged separately"),
//                new ElementRule(UNCONDITIONAL_S_COLUMN, new XMLSyntaxRule[]{
//                       new StringAttributeRule(Attribute.NAME, "The column name")}),
//                new ElementRule(UNCONDITIONAL_N_COLUMN, new XMLSyntaxRule[]{
//                        new StringAttributeRule(Attribute.NAME, "The column name")}),
        };
    };

    final private TraceList traceListN;
    final private TraceList traceListS;
    final private int numSites;
    private OutputFormat format;


    private int fieldWidth;
    private int firstField;
    private NumberFormatter numberFormatter;


}
