/*
 * DnDsPerSiteAnalysis.java
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
import dr.stats.DiscreteStatistics;
import dr.util.*;
import dr.xml.*;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.*;

/**
 * @author Philippe Lemey
 * @author Marc A. Suchard
 */
public class DnDsPerSiteAnalysis implements Citable {
    public static final String DNDS_PER_SITE_ANALYSIS = "dNdSPerSiteAnalysis";
    public static final String BURN_IN = "burnin";
    public static final String CUTOFF = "cutoff";
    public static final String PROPORTION = "proportion";
    public static final String INCLUDE_SIGNIFICANT_SYMBOL = "includeSymbol";
    public static final String INCLUDE_SIGNIFICANCE_LEVEL = "includeLevel";
    public static final String INCLUDE_SITE_CLASSIFICATION = "includeClassification";
    public static final String SIGNIFICANCE_TEST = "test";
    public static final String SEPARATOR_STRING = "separator";
    public static final String INCLUDE_SIMULATION_OUTCOME = "simulationOutcome";
    public static final String INCLUDE_HPD = "includeHPD";
    public static final String INCLUDE_CPD = "includeCPD";
    public static final String SITE_SIMULATION = "siteSimulation";

    public DnDsPerSiteAnalysis(TraceList traceList) {
        this.traceList = traceList;
        this.numSites = traceList.getTraceCount();
        this.format = new OutputFormat();


        fieldWidth = 14;
        firstField = 10;
        numberFormatter = new NumberFormatter(6);
        numberFormatter.setPadding(true);
        numberFormatter.setFieldWidth(fieldWidth);
    }

    public void setIncludeMean(boolean b) {
        format.includeMean = b;
    }

    public void setIncludeHPD(boolean b) {
        format.includeHPD = b;
    }

    public void setIncludeCPD(boolean b) {
        format.includeCPD = b;
    }

    public void setIncludeSignificanceLevel(boolean b) {
        format.includeSignificanceLevel = b;
    }

    public void setIncludeSignificantSymbol(boolean b) {
        format.includeSignificantSymbol = b;
    }

    public void setIncludeSimulationOutcome(boolean b) {
        format.includeSimulationOutcome = b;
    }

    public boolean getIncludeSimulationOutcome() {
        return(format.includeSimulationOutcome);
    }

    public void setProportion(double d) {
        format.proportion = d;
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

    public void setSignificanceTest(SignificanceTest t) {
        format.test = t;
    }

    private String toStringSite(int index, OutputFormat format) {
        StringBuilder sb = new StringBuilder();
        traceList.analyseTrace(index);
        TraceDistribution distribution = traceList.getDistributionStatistics(index);
        sb.append(numberFormatter.formatToFieldWidth(Integer.toString(index + 1), firstField));

        double[] hpd = new double[2];
        double[] cpd = new double[2];

        if (format.proportion == 0.95){
            hpd[0] = distribution.getLowerHPD();
            hpd[1] = distribution.getUpperHPD();
            cpd[0] = distribution.getLowerCPD();
            cpd[1] = distribution.getUpperCPD();
            // this is weird, yes. But we used these 'HPD's/CPDs to obtain the ROC curves for different cut-offs
        }  else if (format.proportion >= 1.0) {
            hpd[0] = cpd[0] = distribution.getMinimum() - (distribution.getMinimum()*(format.proportion - 1.0));
            hpd[1] = cpd[1] = distribution.getMaximum() + (distribution.getMaximum()*(format.proportion - 1.0));
        } else {
//                distribution does not allow to specify proportion
            hpd = getCustomHPDInterval(format.proportion,traceList.getValues(index));
            double[] proportions = new double[]{((1.0-format.proportion)/2.0),(format.proportion+((1.0-format.proportion)/2.0))};
            cpd = getCustomCPDInterval(proportions,traceList.getValues(index));

//                    System.out.println("HPD proportion is "+format.proportion);
        }

        if (format.includeMean) {
            sb.append(format.separator);
            sb.append(numberFormatter.format(distribution.getMean()));
        }
        if (format.includeHPD) {
            sb.append(format.separator);
            sb.append(numberFormatter.format(hpd[0]));
            sb.append(format.separator);
            sb.append(numberFormatter.format(hpd[1]));
        }
        if (format.includeCPD) {
            sb.append(format.separator);
            sb.append(numberFormatter.format(cpd[0]));
            sb.append(format.separator);
            sb.append(numberFormatter.format(cpd[1]));
        }
        if (format.includeSignificanceLevel || format.includeSignificantSymbol || format.includeSiteClassification || format.includeSimulationOutcome) {
            boolean isSignificant = false;
            String classification = "0";
            String level;
            if (format.test == SignificanceTest.NOT_EQUAL_HPD) {
                if (hpd[0] < format.cutoff && hpd[1] < format.cutoff) {
                    level = numberFormatter.formatToFieldWidth(">"+format.proportion, fieldWidth);
                    isSignificant = true;
                    classification = "-";
                } else if (hpd[0] > format.cutoff && hpd[1] > format.cutoff) {
                    level = numberFormatter.formatToFieldWidth(">"+format.proportion, fieldWidth);
                    isSignificant = true;
                    classification = "+";
                } else {
                    level = numberFormatter.formatToFieldWidth("<="+format.proportion, fieldWidth);
                }
            } else if(format.test == SignificanceTest.NOT_EQUAL_CPD) {
                if (cpd[0] < format.cutoff && cpd[1] < format.cutoff) {
                    level = numberFormatter.formatToFieldWidth(">"+format.proportion, fieldWidth);
                    isSignificant = true;
                    classification = "-";
                } else if (cpd[0] > format.cutoff && cpd[1] > format.cutoff) {
                    level = numberFormatter.formatToFieldWidth(">"+format.proportion, fieldWidth);
                    isSignificant = true;
                    classification = "+";
                } else {
                    level = numberFormatter.formatToFieldWidth("<="+format.proportion, fieldWidth);
                }
            } else {
                List values = traceList.getValues(index);
                double levelPosValue = 0.0;
                double levelNegValue = 0.0;
                int total = 0;
                for (Object obj : values) {
                    double d = ((Number) obj).doubleValue();
//                    if ((format.test == SignificanceTest.LESS_THAN && d < format.cutoff) ||
//                            (format.test == SignificanceTest.GREATER_THAN && d > format.cutoff)) {
                    if (d < format.cutoff) {
                        if(format.test == SignificanceTest.LESS_THAN || format.test == SignificanceTest.LESS_OR_GREATER_THAN) {
                            levelNegValue++;
                        }
                    } else if (d > format.cutoff){
                        if (format.test == SignificanceTest.GREATER_THAN || format.test == SignificanceTest.LESS_OR_GREATER_THAN){
                            levelPosValue++;
                        }
                    }
                    total++;
                }
                levelPosValue /= total;
                levelNegValue /= total;
                if (levelPosValue > format.proportion) {
                    isSignificant = true;
                    classification = "+";
                } else if (levelNegValue > format.proportion) {
                    isSignificant = true;
                    classification = "-";
                }
                if (levelPosValue > levelNegValue) {
                    level = numberFormatter.format(levelPosValue);
                                    } else {
                    level = numberFormatter.format(levelNegValue);
                                    }
            }

            if (format.includeSignificanceLevel) {
                sb.append(format.separator);
                sb.append(level);
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
        }
        sb.append("\n");
        return sb.toString();
    }

    public String header(OutputFormat format) {
        StringBuilder sb = new StringBuilder();
        sb.append("# Some information here\n");
        sb.append("# Please cite: " + Citable.Utils.getCitationString(this));


        sb.append(numberFormatter.formatToFieldWidth("Site", firstField));

        if (format.includeMean) {
            sb.append(format.separator);
            sb.append(numberFormatter.formatToFieldWidth("Mean", fieldWidth));
        }

        if (format.includeHPD) {
            sb.append(format.separator);
            sb.append(numberFormatter.formatToFieldWidth("HPD Low", fieldWidth));
            sb.append(format.separator);
            sb.append(numberFormatter.formatToFieldWidth("HPD Up", fieldWidth));
        }

        if (format.includeCPD) {
            sb.append(format.separator);
            sb.append(numberFormatter.formatToFieldWidth("CPD Low", fieldWidth));
            sb.append(format.separator);
            sb.append(numberFormatter.formatToFieldWidth("CPD Up", fieldWidth));
        }

        if (format.includeSignificanceLevel) {
            sb.append(format.separator);
            sb.append(numberFormatter.formatToFieldWidth("Level", fieldWidth));
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
        boolean includeMean;
        boolean includeHPD;
        boolean includeCPD;
        boolean includeSignificanceLevel;
        boolean includeSignificantSymbol;
        boolean includeSiteClassification;
        boolean includeSimulationOutcome;
        String[] siteSimulation;
        double cutoff;
        double proportion;
        SignificanceTest test;
        String separator;

        OutputFormat() {
            this(true, false, true, true, true, true, false, null, 1.0, 0.95, SignificanceTest.NOT_EQUAL_CPD, "\t");
        }

        OutputFormat(boolean includeMean,
                     boolean includeHPD,
                     boolean includeCPD,
                     boolean includeSignificanceLevel,
                     boolean includeSignificantSymbol,
                     boolean includeSiteClassification,
                     boolean includeSimulationOutcome,
                     String[] siteSimulation,
                     double cutoff,
                     double proportion,
                     SignificanceTest test,
                     String separator) {
            this.includeMean = includeMean;
            this.includeHPD = includeHPD;
            this.includeCPD = includeCPD;
            this.includeSignificanceLevel = includeSignificanceLevel;
            this.includeSignificantSymbol = includeSignificantSymbol;
            this.includeSiteClassification = includeSiteClassification;
            this.includeSimulationOutcome = includeSimulationOutcome;
            this.siteSimulation = siteSimulation;
            this.cutoff = cutoff;
            this.proportion = proportion;
            this.test = test;
            this.separator = separator;
        }
    }

    public enum SignificanceTest {
        GREATER_THAN("gt"),    //>
        LESS_THAN("lt"),       //<
        NOT_EQUAL_HPD("ne_HPD"),       //!=
        NOT_EQUAL_CPD("ne_CPD"),       //!=
        LESS_OR_GREATER_THAN("logt"); //<>

        private SignificanceTest(String text) {
            this.text = text;
        }

        public String getText() {
            return text;
        }

        public static SignificanceTest parseFromString(String text) {
            for (SignificanceTest test : SignificanceTest.values()) {
                if (test.getText().compareToIgnoreCase(text) == 0)
                    return test;
            }
            return null;
        }

        private final String text;
    }

    private static double[] getCustomHPDInterval(double proportion, List list) {

        double returnArray[] = new double[2];
        int length = list.size();
        int[] indices = new int[length];
        Double[] resultObjArray = (Double[]) list.toArray( new Double[0] );
        double[] result = toPrimitiveDoubleArray(resultObjArray);
        HeapSort.sort(result, indices);
        double minRange = Double.MAX_VALUE;
        int hpdIndex = 0;

        int diff = (int)Math.round(proportion * (double)length);

        for (int i = 0; i <= (length - diff); i++) {
            double minValue = result[indices[i]];
            double maxValue = result[indices[i+diff-1]];
            double range = Math.abs(maxValue - minValue);
            if (range < minRange) {
                minRange = range;
                hpdIndex = i;
            }
        }
        returnArray[0] = result[indices[hpdIndex]];
        returnArray[1] = result[indices[hpdIndex+diff-1]];
        return returnArray;
    }

    private static double[] getCustomCPDInterval(double[] proportions, List list) {

        double returnArray[] = new double[2];
        int length = list.size();
        int[] indices = new int[length];
        Double[] resultObjArray = (Double[]) list.toArray( new Double[0] );
        double[] result = toPrimitiveDoubleArray(resultObjArray);
        HeapSort.sort(result, indices);

        returnArray[0] = DiscreteStatistics.quantile(proportions[0], result, indices);
        returnArray[1] = DiscreteStatistics.quantile(proportions[1], result, indices);

        return returnArray;
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
            return DNDS_PER_SITE_ANALYSIS;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            String fileName = xo.getStringAttribute(FileHelpers.FILE_NAME);
            try {

                File file = new File(fileName);
                String name = file.getName();
                String parent = file.getParent();

                if (!file.isAbsolute()) {
                    parent = System.getProperty("user.dir");
                }

                file = new File(parent, name);

                fileName = file.getAbsolutePath();

                LogFileTraces traces = new LogFileTraces(fileName, file);
                traces.loadTraces();
                long maxState = traces.getMaxState();

                // leaving the burnin attribute off will result in 10% being used
                long burnin = xo.getAttribute(BURN_IN, maxState / 10);
                //TODO: implement custom burn-in

                if (burnin < 0 || burnin >= maxState) {
                    burnin = maxState / 5;
                    System.out.println("WARNING: Burn-in larger than total number of states - using 20%");
                }

                traces.setBurnIn(burnin);

                // TODO: Filter traces to include only dNdS columns

                DnDsPerSiteAnalysis analysis = new DnDsPerSiteAnalysis(traces);

                analysis.setSignificanceTest(
                        SignificanceTest.parseFromString(
                                xo.getAttribute(SIGNIFICANCE_TEST, SignificanceTest.NOT_EQUAL_CPD.getText())
                        )
                );
                analysis.setCutoff(xo.getAttribute(CUTOFF, 1.0));
                analysis.setProportion(xo.getAttribute(PROPORTION, 0.95));
                analysis.setSeparator(xo.getAttribute(SEPARATOR_STRING, "\t"));
                analysis.setIncludeHPD(xo.getAttribute(INCLUDE_HPD, false));
                analysis.setIncludeCPD(xo.getAttribute(INCLUDE_HPD, true));
                analysis.setIncludeSignificanceLevel(xo.getAttribute(INCLUDE_SIGNIFICANCE_LEVEL, false));
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
            return "Performs a trace dN/dS analysis.";
        }

        public Class getReturnType() {
            return DnDsPerSiteAnalysis.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                AttributeRule.newDoubleRule(CUTOFF, true),
                AttributeRule.newDoubleRule(PROPORTION, true),
                AttributeRule.newIntegerRule(BURN_IN,true),
                AttributeRule.newBooleanRule(INCLUDE_HPD, true),
                AttributeRule.newBooleanRule(INCLUDE_CPD, true),
                AttributeRule.newBooleanRule(INCLUDE_SIGNIFICANT_SYMBOL, true),
                AttributeRule.newBooleanRule(INCLUDE_SIGNIFICANCE_LEVEL, true),
                AttributeRule.newBooleanRule(INCLUDE_SITE_CLASSIFICATION, true),
                AttributeRule.newBooleanRule(INCLUDE_SIMULATION_OUTCOME, true),
                AttributeRule.newStringRule(SITE_SIMULATION, true),
                AttributeRule.newStringRule(SIGNIFICANCE_TEST, true),
                AttributeRule.newStringRule(SEPARATOR_STRING, true),
                new StringAttributeRule(FileHelpers.FILE_NAME,
                        "The traceName of a BEAST log file (can not include trees, which should be logged separately"),
//                new ElementRule(UNCONDITIONAL_S_COLUMN, new XMLSyntaxRule[]{
//                       new StringAttributeRule(Attribute.NAME, "The column name")}),
//                new ElementRule(UNCONDITIONAL_N_COLUMN, new XMLSyntaxRule[]{
//                        new StringAttributeRule(Attribute.NAME, "The column name")}),
        };
    };

    final private TraceList traceList;
    final private int numSites;
    private OutputFormat format;


    private int fieldWidth;
    private int firstField;
    private NumberFormatter numberFormatter;


//    private String separator = "\t";
//    final private int numSamples;

//    private double[][][] allSamples;
//    final private static int NUM_VARIABLES = 4;
//    final private static int COND_S = 0;
//    final private static int UNCOND_S = 1;
//    final private static int COND_N = 2;
//    final private static int UNCOND_N = 3;
//    final private static String[] names = {COND_SPERSITE_COLUMNS, UNCOND_SPERSITE_COLUMNS, COND_NPERSITE_COLUMNS, UNCOND_NPERSITE_COLUMNS};
//    private double[][][] smoothSamples;
//    private double[][][] smoothDnDsSamples;

    private static final boolean DEBUG = true;

//    private double[][] rawMeanStats;
//    private double[][] smoothMeanStats;
//    private double[][] smoothMeanDnDsStats;
//    private double[][][] smoothHPDDnDsStats;

}
