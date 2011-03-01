package dr.inference.trace;

import dr.inferencexml.trace.MarginalLikelihoodAnalysisParser;
import dr.util.*;
import dr.xml.*;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

/**
 * @author Philippe Lemey
 * @author Marc A. Suchard
 */
public class DnDsPerSiteAnalysis implements Citable {
    public static final String DNDS_PER_SITE_ANALYSIS = "dNdSPerSiteAnalysis";
    public static final String CUTOFF = "cutoff";
    public static final String PROPORTION = "proportion";
    public static final String INCLUDE_SIGNIFICANT_SYMBOL = "includeSymbol";
    public static final String INCLUDE_SIGNIFICANCE_LEVEL = "includeLevel";
    public static final String INCLUDE_SITE_CLASSIFICATION = "includeClassification";
    public static final String SIGNIFICANCE_TEST = "test";
    public static final String SEPARATOR_STRING = "separator";
    public static final String INCLUDE_SIMULATION_OUTCOME = "simulationOutcome";
    public static final String INCLUDE_HPD = "includeHPD";

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

    public void setIncludeSignificanceLevel(boolean b) {
        format.includeSignificanceLevel = b;
    }

    public void setIncludeSignificantSymbol(boolean b) {
        format.includeSignificantSymbol = b;
    }

    public void setIncludeSimulationOutcome(boolean b) {
        format.includeSimulationOutcome = b;
    }

    public void setProportion(double d) {
        format.proportion = d;
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
        if (format.includeMean) {
            sb.append(format.separator);
            sb.append(numberFormatter.format(distribution.getMean()));
        }
        if (format.includeHPD) {
            sb.append(format.separator);
            sb.append(numberFormatter.format(distribution.getLowerHPD()));
            sb.append(format.separator);
            sb.append(numberFormatter.format(distribution.getUpperHPD()));
        }
        if (format.includeSignificanceLevel || format.includeSignificantSymbol || format.includeSiteClassification || format.includeSimulationOutcome) {
            boolean isSignificant = false;
            String classification = "0";
            String level;
            if (format.test == SignificanceTest.NOT_EQUAL) {
                double[] hpd = getHPDInterval(format.proportion,traceList.getValues(index));
//                distribution does not allow to specify proportion
//                double lower = distribution.getLowerHPD();
//                double upper = distribution.getUpperHPD();

//                if ((lower < format.cutoff && upper < format.cutoff) ||
//                        (lower > format.cutoff && upper > format.cutoff)) {
                if (hpd[0] < format.cutoff && hpd[1] < format.cutoff) {
                    level = numberFormatter.formatToFieldWidth(">0.95", fieldWidth);
                    isSignificant = true;
                    classification = "-";
                } else if (hpd[0] > format.cutoff && hpd[1] > format.cutoff) {
                    level = numberFormatter.formatToFieldWidth(">0.95", fieldWidth);
                    isSignificant = true;
                    classification = "+";
                } else {
                    level = numberFormatter.formatToFieldWidth("<=0.95", fieldWidth);
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
                sb.append(simulated[index]);
                sb.append(format.separator);
                if (simulated[index].equals("+") || simulated[index].equals("-")) {
                    if (classification.equals(simulated[index])){
                        sb.append("TP");   // True Positive
                    } else {
                        sb.append("FN");   // True Negative
                    }
                }  else {
                    if (classification.equals(simulated[index])){
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

        return sb.toString();
    }

    public List<Citation> getCitations() {
        List<Citation> citations = new ArrayList<Citation>();
        citations.add(
                new Citation(
                        new Author[]{
                                new Author("P", "Lemey"),
                                new Author("VN", "Minin"),
                                new Author("MA", "Suchard")
                        },
                        Citation.Status.IN_PREPARATION
                )
        );
        return citations;
    }

    private class OutputFormat {
        boolean includeMean;
        boolean includeHPD;
        boolean includeSignificanceLevel;
        boolean includeSignificantSymbol;
        boolean includeSiteClassification;
        boolean includeSimulationOutcome;
        double cutoff;
        double proportion;
        SignificanceTest test;
        String separator;

        OutputFormat() {
            this(true, true, true, true, true, false, 1.0, 0.95, SignificanceTest.NOT_EQUAL, "\t");
        }

        OutputFormat(boolean includeMean,
                     boolean includeHPD,
                     boolean includeSignificanceLevel,
                     boolean includeSignificantSymbol,
                     boolean includeSiteClassification,
                     boolean includeSimulationOutcome,
                     double cutoff,
                     double proportion,
                     SignificanceTest test,
                     String separator) {
            this.includeMean = includeMean;
            this.includeHPD = includeHPD;
            this.includeSignificanceLevel = includeSignificanceLevel;
            this.includeSignificantSymbol = includeSignificantSymbol;
            this.includeSiteClassification = includeSiteClassification;
            this.includeSimulationOutcome = includeSimulationOutcome;
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

    private static double[] getHPDInterval(double proportion, List list) {

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

    private static double[] toPrimitiveDoubleArray(Double[] array){
        double[] returnArray = new double[array.length];
        for(int i = 0; i < array.length; i++ ){
            returnArray[i] = array[i].doubleValue();
        }
        return returnArray;
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
                int maxState = traces.getMaxState();

                // leaving the burnin attribute off will result in 10% being used
                int burnin = xo.getAttribute(MarginalLikelihoodAnalysisParser.BURN_IN, maxState / 10);
                //TODO: implement custom burn-in

                if (burnin < 0 || burnin >= maxState) {
                    burnin = maxState / 5;
                    System.out.println("WARNING: Burn-in larger than total number of states - using to 20%");
                }

                traces.setBurnIn(burnin);

                // TODO: Filter traces to include only dNdS columns

                DnDsPerSiteAnalysis analysis = new DnDsPerSiteAnalysis(traces);

                analysis.setSignificanceTest(
                        SignificanceTest.parseFromString(
                                xo.getAttribute(SIGNIFICANCE_TEST, SignificanceTest.NOT_EQUAL.getText())
                        )
                );
                analysis.setCutoff(xo.getAttribute(CUTOFF, 1.0));
                analysis.setProportion(xo.getAttribute(PROPORTION, 0.95));
                analysis.setSeparator(xo.getAttribute(SEPARATOR_STRING, "\t"));
                analysis.setIncludeHPD(xo.getAttribute(INCLUDE_HPD, true));
                analysis.setIncludeSignificanceLevel(xo.getAttribute(INCLUDE_SIGNIFICANCE_LEVEL, false));
                analysis.setIncludeSignificantSymbol(xo.getAttribute(INCLUDE_SIGNIFICANT_SYMBOL, true));
                analysis.setIncludeSiteClassification(xo.getAttribute(INCLUDE_SITE_CLASSIFICATION, true));
                analysis.setIncludeSimulationOutcome(xo.getAttribute(INCLUDE_SIMULATION_OUTCOME, false));

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
                AttributeRule.newDoubleArrayRule(CUTOFF, true),
                AttributeRule.newDoubleArrayRule(PROPORTION, true),
                AttributeRule.newBooleanRule(INCLUDE_HPD, true),
                AttributeRule.newBooleanRule(INCLUDE_SIGNIFICANT_SYMBOL, true),
                AttributeRule.newBooleanRule(INCLUDE_SIGNIFICANCE_LEVEL, true),
                AttributeRule.newBooleanRule(INCLUDE_SITE_CLASSIFICATION, true),
                AttributeRule.newBooleanRule(INCLUDE_SIMULATION_OUTCOME, true),
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
    private String simulated[] = {"-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "+", "+", "+", "+", "+", "+", "+", "+", "+", "+", "+", "+", "+", "+", "+", "+", "+", "+", "+", "+", "+", "+", "+", "+", "+", "+", "+", "+", "+", "+", "+", "+", "+", "+", "+", "+", "+", "+", "+", "+", "+", "+", "+", "+", "+", "+", "+", "+", "+", "+", "+", "+", "+", "+", "+", "+", "+", "+", "+", "+", "+", "+", "+", "+", "+", "+", "+", "+", "+", "+", "+", "+", "+", "+", "+"};


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
