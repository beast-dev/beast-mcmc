/*
 * CnCsToDnDsPerSiteAnalysis.java
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

import dr.math.EmpiricalBayesPoissonSmoother;
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
public class CnCsToDnDsPerSiteAnalysis implements Citable {

    public static final String CNCS_TO_DNDS_PER_SITE_ANALYSIS = "cNcSTodNdSPerSiteAnalysis";
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
    public static final String SITE_SIMULATION = "siteSimulation";
    public static final String CN = "CN";
    public static final String CS = "CS";
    public static final String USE_SAMPLE = "sample";

    public CnCsToDnDsPerSiteAnalysis(TraceList traceListN, TraceList traceListS) {
        this.traceListN = traceListN;
        this.traceListS = traceListS;
        this.numSites = (traceListN.getTraceCount()) / 2;
        this.format = new OutputFormat();


        setUseSample(false);


        double[][] allCn = new double[numSites][this.traceListN.getStateCount()];
        double[][] allUn = new double[numSites][this.traceListN.getStateCount()];
        double[][] allCs = new double[numSites][this.traceListS.getStateCount()];
        double[][] allUs = new double[numSites][this.traceListS.getStateCount()];

        for (int o = 0; o < numSites; o++) {
            allCn[o] = listToDoubleArray(traceListN.getValues(o));
            allUn[o] = listToDoubleArray(traceListN.getValues(numSites + o));
            allCs[o] = listToDoubleArray(traceListS.getValues(o));
            allUs[o] = listToDoubleArray(traceListS.getValues(numSites + o));
        }

        allCn = transpose(allCn);
        allUn = transpose(allUn);
        allCs = transpose(allCs);
        allUs = transpose(allUs);

        double[][] tempAllSmoothedCn = new double[this.traceListN.getStateCount()][numSites];
        double[][] tempAllSmoothedUn = new double[this.traceListN.getStateCount()][numSites];
        double[][] tempAllSmoothedCs = new double[this.traceListS.getStateCount()][numSites];
        double[][] tempAllSmoothedUs = new double[this.traceListS.getStateCount()][numSites];

        boolean first = true;

        for (int p = 0; p < this.traceListN.getStateCount(); p++) {

            // Testing code:

//            if (first) {
//                tempAllSmoothedCn[p] = EmpiricalBayesPoissonSmoother.smooth(allCn[p]);
//                System.err.println("Smooth values: " + new Vector(tempAllSmoothedCn[p]));
//
//                tempAllSmoothedCn[p] = EmpiricalBayesPoissonSmoother.smoothWithSample(allCn[p]);
//                System.err.println("Sample values: " + new Vector(tempAllSmoothedCn[p]));
//                first = false;
//            }

//            System.exit(-1);


            if (format.useSample) {
                tempAllSmoothedCn[p] = EmpiricalBayesPoissonSmoother.smoothWithSample(allCn[p]);
                tempAllSmoothedUn[p] = EmpiricalBayesPoissonSmoother.smoothWithSample(allUn[p]);
                tempAllSmoothedCs[p] = EmpiricalBayesPoissonSmoother.smoothWithSample(allCs[p]);
                tempAllSmoothedUs[p] = EmpiricalBayesPoissonSmoother.smoothWithSample(allUs[p]);
            } else {
                tempAllSmoothedCn[p] = EmpiricalBayesPoissonSmoother.smooth(allCn[p]);
                tempAllSmoothedUn[p] = EmpiricalBayesPoissonSmoother.smooth(allUn[p]);
                tempAllSmoothedCs[p] = EmpiricalBayesPoissonSmoother.smooth(allCs[p]);
                tempAllSmoothedUs[p] = EmpiricalBayesPoissonSmoother.smooth(allUs[p]);
            }

        }

        allSmoothedCn = transpose(tempAllSmoothedCn);
        allSmoothedUn = transpose(tempAllSmoothedUn);
        allSmoothedCs = transpose(tempAllSmoothedCs);
        allSmoothedUs = transpose(tempAllSmoothedUs);

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
        return (format.includeSimulationOutcome);
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

    public void setUseSample(boolean u) {
        format.useSample = u;
    }

    private String toStringSite(int index, OutputFormat format) {

        double[] dN = getRatioArray(allSmoothedCn[index], allSmoothedUn[index]);
        double[] dS = getRatioArray(allSmoothedCs[index], allSmoothedUs[index]);

        double[] omegas = getRatioArray(dN, dS);
//        for (int x = 0; x < dN.length; x ++){
//            System.out.println(index+"\t"+allSmoothedCn[index][x]+"\t"+allSmoothedUn[index][x]+"\t"+dN[x]+"\t"+allSmoothedCs[index][x]+"\t"+allSmoothedUs[index][x]+"\t"+dS[x]+"\t"+omegas[x]);
//        }

        StringBuilder sb = new StringBuilder();
        sb.append(numberFormatter.formatToFieldWidth(Integer.toString(index + 1), firstField));

        double[] hpd = new double[2];
        double[] minMax = getMinMax(omegas);
        // this is weird, yes. But we used these 'HPD's to obtain the ROC curves for different cut-offs
        if (format.proportion >= 1.0) {
            hpd[0] = minMax[0] - (minMax[0] * (format.proportion - 1.0));
            hpd[1] = minMax[1] + (minMax[1] * (format.proportion - 1.0));
            System.out.println("hpd = " + hpd[0] + " - " + hpd[1]);
        } else {
            hpd = getHPDInterval(format.proportion, omegas);
        }

        if (format.includeMean) {
            sb.append(format.separator);
//            sb.append(numberFormatter.format(DiscreteStatistics.mean(omegas)));
            sb.append(numberFormatter.format(DiscreteStatistics.median(omegas)));
        }
        if (format.includeHPD) {
            sb.append(format.separator);
            sb.append(numberFormatter.format(hpd[0]));
            sb.append(format.separator);
            sb.append(numberFormatter.format(hpd[1]));
        }
        if (format.includeSignificanceLevel || format.includeSignificantSymbol || format.includeSiteClassification || format.includeSimulationOutcome) {
            boolean isSignificant = false;
            String classification = "0";
            String level;
            if (format.test == SignificanceTest.NOT_EQUAL) {
                if (hpd[0] < format.cutoff && hpd[1] < format.cutoff) {
                    level = numberFormatter.formatToFieldWidth(">" + format.proportion, fieldWidth);
                    isSignificant = true;
                    classification = "-";
                } else if (hpd[0] > format.cutoff && hpd[1] > format.cutoff) {
                    level = numberFormatter.formatToFieldWidth(">" + format.proportion, fieldWidth);
                    isSignificant = true;
                    classification = "+";
                } else {
                    level = numberFormatter.formatToFieldWidth("<=" + format.proportion, fieldWidth);
                }
            } else {
                double levelPosValue = 0.0;
                double levelNegValue = 0.0;
                int total = 0;
                for (double w : omegas) {
//                    if ((format.test == SignificanceTest.LESS_THAN && d < format.cutoff) ||
//                            (format.test == SignificanceTest.GREATER_THAN && d > format.cutoff)) {
                    if (w < format.cutoff) {
                        if (format.test == SignificanceTest.LESS_THAN || format.test == SignificanceTest.LESS_OR_GREATER_THAN) {
                            levelNegValue++;
                        }
                    } else if (w > format.cutoff) {
                        if (format.test == SignificanceTest.GREATER_THAN || format.test == SignificanceTest.LESS_OR_GREATER_THAN) {
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
                    if (classification.equals(format.siteSimulation[index])) {
                        sb.append("TP");   // True Positive
                    } else {
                        sb.append("FN");   // True Negative
                    }
                } else {
                    if (classification.equals(format.siteSimulation[index])) {
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
        sb.append("# Please cite: " + Utils.getCitationString(this));


        sb.append(numberFormatter.formatToFieldWidth("Site", firstField));

        if (format.includeMean) {
            sb.append(format.separator);
            sb.append(numberFormatter.formatToFieldWidth("Mean", fieldWidth));
        }

        if (format.includeHPD) {
            sb.append(format.separator);
            sb.append(numberFormatter.formatToFieldWidth("Lower", fieldWidth));
            sb.append(format.separator);
            sb.append(numberFormatter.formatToFieldWidth("Upper", fieldWidth));
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

        sb.append(toStringSite(0, format));

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
        boolean useSample;
        boolean includeMean;
        boolean includeHPD;
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
            this(false, true, true, true, true, true, false, null, 1.0, 0.95, SignificanceTest.NOT_EQUAL, "\t");
        }

        OutputFormat(boolean useSample,
                     boolean includeMean,
                     boolean includeHPD,
                     boolean includeSignificanceLevel,
                     boolean includeSignificantSymbol,
                     boolean includeSiteClassification,
                     boolean includeSimulationOutcome,
                     String[] siteSimulation,
                     double cutoff,
                     double proportion,
                     SignificanceTest test,
                     String separator) {
            this.useSample = useSample;
            this.includeMean = includeMean;
            this.includeHPD = includeHPD;
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
        NOT_EQUAL("ne"),       //!=
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

    private static double[] listToDoubleArray(List list) {
        Double[] resultObjArray = (Double[]) list.toArray(new Double[0]);
        double[] result = toPrimitiveDoubleArray(resultObjArray);
        return result;
    }

    private double[][] transpose(double[][] in) {
        double[][] out = new double[in[0].length][in.length];
        for (int r = 0; r < in.length; r++) {
            for (int c = 0; c < in[0].length; c++) {
                out[c][r] = in[r][c];
            }
        }
        return out;
    }

    private static double[] getFirstHalfArray(double[] condAndUncond) {
        int count = (condAndUncond.length) / 2;
        double[] returnArray = new double[count];
        for (int a = 0; a < count; a++) {
            returnArray[a] = condAndUncond[a];
        }
        return returnArray;
    }

    private static double[] getSecondHalfArray(double[] condAndUncond) {
        int count = (condAndUncond.length) / 2;
        double[] returnArray = new double[count];
        for (int a = 0; a < count; a++) {
            returnArray[a] = condAndUncond[count + a];
        }
        return returnArray;
    }

    private static double[] getRatioArray(double[] enumArray, double[] denomArray) {
        double[] returnArray = new double[enumArray.length];
        for (int x = 0; x < enumArray.length; x++) {
            returnArray[x] = enumArray[x] / denomArray[x];
        }
        return returnArray;
    }

    private static double[] getMinMax(double[] values) {
        double[] returnArray = new double[2];

        double min = Double.MAX_VALUE;
        double max = 0;

        for (int x = 0; x < values.length; x++) {
            if (values[x] > max) {
                max = values[x];
            }
            if (values[x] < min) {
                min = values[x];
            }
        }

        returnArray[0] = min;
        returnArray[1] = max;

        return returnArray;

    }

    private static double[] getHPDInterval(double proportion, double[] values) {

        double[] returnArray = new double[2];

        int length = values.length;
        int[] indices = new int[length];
        HeapSort.sort(values, indices);
        double minRange = Double.MAX_VALUE;
        int hpdIndex = 0;

        int diff = (int) Math.round(proportion * (double) length);

        for (int i = 0; i <= (length - diff); i++) {
            double minValue = values[indices[i]];
            double maxValue = values[indices[i + diff - 1]];
            double range = Math.abs(maxValue - minValue);
            if (range < minRange) {
                minRange = range;
                hpdIndex = i;
            }
        }
        returnArray[0] = values[indices[hpdIndex]];
        returnArray[1] = values[indices[hpdIndex + diff - 1]];
        return returnArray;
    }

    private static double[] toPrimitiveDoubleArray(Double[] array) {
        double[] returnArray = new double[array.length];
        for (int i = 0; i < array.length; i++) {
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
            return CNCS_TO_DNDS_PER_SITE_ANALYSIS;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            String fileNameCN = xo.getStringAttribute(FileHelpers.FILE_NAME + CN);
            String fileNameCS = xo.getStringAttribute(FileHelpers.FILE_NAME + CS);
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

                fileNameCN = fileCN.getAbsolutePath();
                fileNameCS = fileCS.getAbsolutePath();

                LogFileTraces tracesCN = new LogFileTraces(fileNameCN, fileCN);
                LogFileTraces tracesCS = new LogFileTraces(fileNameCS, fileCS);
                tracesCN.loadTraces();
                tracesCS.loadTraces();

                long maxStateCN = tracesCN.getMaxState();
                long maxStateCS = tracesCS.getMaxState();

                if (maxStateCN != maxStateCS) {
                    System.err.println("max states in" + fileNameCN + " and " + fileNameCS + " are not equal");
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

                CnCsToDnDsPerSiteAnalysis analysis = new CnCsToDnDsPerSiteAnalysis(tracesCN, tracesCS);

                analysis.setCutoff(xo.getAttribute(CUTOFF, 1.0));
                analysis.setProportion(xo.getAttribute(PROPORTION, 0.95));
                analysis.setSeparator(xo.getAttribute(SEPARATOR_STRING, "\t"));
                analysis.setUseSample(xo.getAttribute(USE_SAMPLE, false));
                analysis.setIncludeHPD(xo.getAttribute(INCLUDE_HPD, true));
                analysis.setIncludeSignificanceLevel(xo.getAttribute(INCLUDE_SIGNIFICANCE_LEVEL, false));
                analysis.setIncludeSignificantSymbol(xo.getAttribute(INCLUDE_SIGNIFICANT_SYMBOL, true));
                analysis.setIncludeSiteClassification(xo.getAttribute(INCLUDE_SITE_CLASSIFICATION, true));
                analysis.setIncludeSimulationOutcome(xo.getAttribute(INCLUDE_SIMULATION_OUTCOME, false));
                if (analysis.getIncludeSimulationOutcome()) {
                    String sites = (String) xo.getAttribute(SITE_SIMULATION, "empty");
                    if (sites.equals("empty")) {
                        System.err.println("you want simulation evaluation but do not provide a site simulation string??");
                    } else {
                        String[] siteSimulation = parseVariableLengthStringArray(sites);
                        analysis.setSiteSimulation(siteSimulation);
                    }
                }

                return analysis;

            } catch (FileNotFoundException fnfe) {
                throw new XMLParseException("File '" + fileNameCN + "  and " + fileNameCS + "' can not be opened for " + getParserName() + " element.");
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
                AttributeRule.newStringRule(FileHelpers.FILE_NAME + CN, true),
                AttributeRule.newStringRule(FileHelpers.FILE_NAME + CS, true),
                AttributeRule.newDoubleRule(CUTOFF, true),
                AttributeRule.newDoubleRule(PROPORTION, true),
                AttributeRule.newIntegerRule(BURN_IN, true),
                AttributeRule.newBooleanRule(USE_SAMPLE, true),
                AttributeRule.newBooleanRule(INCLUDE_HPD, true),
                AttributeRule.newBooleanRule(INCLUDE_SIGNIFICANT_SYMBOL, true),
                AttributeRule.newBooleanRule(INCLUDE_SIGNIFICANCE_LEVEL, true),
                AttributeRule.newBooleanRule(INCLUDE_SITE_CLASSIFICATION, true),
                AttributeRule.newBooleanRule(INCLUDE_SIMULATION_OUTCOME, true),
                AttributeRule.newStringRule(SITE_SIMULATION, true),
                AttributeRule.newStringRule(SIGNIFICANCE_TEST, true),
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
    final private double[][] allSmoothedCn;
    final private double[][] allSmoothedUn;
    final private double[][] allSmoothedCs;
    final private double[][] allSmoothedUs;
    private OutputFormat format;


    private int fieldWidth;
    private int firstField;
    private NumberFormatter numberFormatter;


}
