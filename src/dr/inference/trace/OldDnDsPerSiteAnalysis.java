/*
 * OldDnDsPerSiteAnalysis.java
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
import dr.math.EmpiricalBayesPoissonSmoother;
import dr.stats.DiscreteStatistics;
import dr.util.Attribute;
import dr.util.FileHelpers;
import dr.util.HeapSort;
import dr.xml.*;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.List;

/**
 * @author Philippe Lemey
 * @author Marc A. Suchard
 */
public class OldDnDsPerSiteAnalysis {
    public static final String DNDS_PERSITE_ANALYSIS = "olddNdSPerSiteAnalysis";
    public static final String COND_SPERSITE_COLUMNS = "conditionalSperSite";
    public static final String UNCOND_SPERSITE_COLUMNS = "unconditionalSperSite";
    public static final String COND_NPERSITE_COLUMNS = "conditionalNperSite";
    public static final String UNCOND_NPERSITE_COLUMNS = "unconditionalNperSite";
    //public static final String ALIGNMENT = "alignment";

    OldDnDsPerSiteAnalysis(double[][] sampleSperSite, double[][] unconditionalS,
                        double[][] sampleNperSite, double[][] unconditionalN) {

        numSites = sampleNperSite.length;
        numSamples = sampleNperSite[0].length;

        allSamples = new double[NUM_VARIABLES][numSites][numSamples];
        allSamples[COND_S] = sampleSperSite;
        allSamples[UNCOND_S] = unconditionalS;
        allSamples[COND_N] = sampleNperSite;
        allSamples[UNCOND_N] = unconditionalN;

        if (DEBUG) {
            System.err.println("sumSites = " + numSites);
            System.err.println("numSamples = " + numSamples);
        }

        // TODO Assert that all arrays have the same dimensions, else throw exception

        smoothSamples = performSmoothing(allSamples);
        // collect 'dN'(=N/U_N)', 'dS', 'dN/dS'
        smoothDnDsSamples = getDnDsSamples(smoothSamples);

        rawMeanStats = computeMeanStats(allSamples);
        smoothMeanStats = computeMeanStats(smoothSamples);
        smoothMeanDnDsStats = computeMeanStats(smoothDnDsSamples);
        smoothHPDDnDsStats = computeHPDStats(smoothDnDsSamples);

    }

//    public double[][] getdSPerSiteSample() {
//        if (!dNAnddSPerSiteSampleCollected) {
//            collectDnAndDsPerSite();
//        }
//        return dSPerSiteSample;
//    }
//
//    public double[][] getdNPerSiteSample() {
//        if (!dNAnddSPerSiteSampleCollected) {
//            collectDnAndDsPerSite();
//        }
//        return dNPerSiteSample;
//    }
//
//    public double[][] getdNdSRatioPerSiteSample() {
//        if (!dNAnddSPerSiteSampleCollected) {
//            collectDnAndDsPerSite();
//        }
//        return dNdSRatioPerSiteSample;
//    }
//
//    public double[][] getEBdNdSratioPerSiteSample() {
//        if (!dNAnddSPerSiteSampleCollected) {
//            collectDnAndDsPerSite();
//        }
//        return EBdNdSratioPerSiteSample;
//    }

    private double[][][] getDnDsSamples(double[][][] smoothedSamples) {
        double[][][] dNdSArray = new double[3][][];
        //get smooth(S)/smooth(u_S)
        dNdSArray[0] = get2DArrayRatio(smoothedSamples[0],smoothedSamples[1]);
        //get smooth(N)/smooth(u_N)
        dNdSArray[1] = get2DArrayRatio(smoothedSamples[2],smoothedSamples[3]);
        //get (smooth(N)/smooth(u_N))/(smooth(S)/smooth(u_S))
        dNdSArray[2] = get2DArrayRatio(dNdSArray[1],dNdSArray[0]);
        return dNdSArray;
    }

    private double[][] get2DArrayRatio (double[][] numerator, double[][] denominator) {
        double[][] returnArray = new double[numerator.length][numerator[0].length];
        for(int site = 0; site < numerator.length; site ++) {
            returnArray[site] = get1DArrayRatio(numerator[site],denominator[site]);
        }
        return returnArray;
    }

    private double[] get1DArrayRatio (double[] numerator, double[] denominator) {
        double[] returnArray = new double[numerator.length];
        for(int sample = 0; sample < numerator.length; sample ++) {
            returnArray[sample] = numerator[sample]/denominator[sample];
        }
        return returnArray;
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

    private double computeDnDs(double[][][] allSamples, int site, int sample) {
        return (allSamples[COND_N][site][sample] / allSamples[COND_S][site][sample]) /
                (allSamples[UNCOND_N][site][sample] / allSamples[UNCOND_S][site][sample]);
    }

    private double[][][] performSmoothing(double[][][] allSamples) {
        double[][][] smoothedArray = new double[allSamples.length][][];
        double[][][] transpose = new double[allSamples.length][][];
        for (int i = 0; i < allSamples.length; i++) {
            transpose[i] = transpose(allSamples[i]); // transpose means: [allSamples][site]
            for (int sample = 0; sample < numSamples; ++sample) {
                transpose[i][sample] = EmpiricalBayesPoissonSmoother.smooth(transpose[i][sample]);
            }
            smoothedArray[i] = transpose(transpose[i]);
        }
        return smoothedArray;

    }

    private double computeRatio(double[][][] allSamples, int site, int sample, int numerator, int denominator) {
        return allSamples[numerator][site][sample] / allSamples[denominator][site][sample];
    }

    private double[][] computeMeanStats(double[][][] allSamples) {
        double[][] statistics = new double[allSamples.length][allSamples[0].length];
        for (int variable = 0; variable < allSamples.length; ++variable) {
            statistics[variable] = mean(allSamples[variable]);
        }
        return statistics;
    }

    private double[][][] computeHPDStats(double[][][] allSamples) {
        // we collect mean and 95% HPD boundaries
        double[][][] statistics = new double[allSamples.length][allSamples[0].length][2];
        for (int variable = 0; variable < allSamples.length; variable++) {
            statistics[variable] = getArrayHPDintervals(allSamples[variable]);
        }
        return statistics;
    }

//    private void collectDnAndDsPerSite() {
//
//        // If output arrays not yet constructed, then allocate memory once
//        if (dSPerSiteSample == null) {
//            dSPerSiteSample = new double[numSites][numSamples];
//            dNPerSiteSample = new double[numSites][numSamples];
//            dNdSRatioPerSiteSample = new double[numSites][numSamples];
//            EBdNdSratioPerSiteSample = new double[numSites][numSamples];
//        }
//
//        // Compute statistics on unsmoothed values
//        for (int site = 0; site < numSites; ++site) {
//            for (int sample = 0; sample < numSamples; ++sample) {
//                dSPerSiteSample[site][sample] = computeRatio(allSamples, site, sample, COND_S, UNCOND_S);
//                dNPerSiteSample[site][sample] = computeRatio(allSamples, site, sample, COND_N, UNCOND_N);
//                dNdSRatioPerSiteSample[site][sample] = computeDnDs(allSamples, site, sample);
//            }
//        }
//
//        if (DEBUG) {
//            System.err.println("allSamples = " + new Vector(allSamples[0][0]));
////            System.exit(-1);
//            System.err.println("dS = " + new Vector(dSPerSiteSample[0]));
//        }
//
//        performSmoothing(allSamples);
//
//        // Compute test statistics on smoothed values
//        for (int site = 0; site < numSites; ++site) {
//            for (int sample = 0; sample < numSamples; ++sample) {
//                dNdSRatioPerSiteSample[site][sample] = computeDnDs(allSamples, site, sample);
//            }
//        }
//
//        dNAnddSPerSiteSampleCollected = true;
//    }


//    public static double[] smooth(double[] in) {
//        final int length = in.length;
//        double[] out = new double[length];
//        double[] gammaStats = getNegBin(in);
//        for (int i = 0; i < length; i++) {
//            out[i] = (in[i] + gammaStats[0]) / (1 + 1 / gammaStats[1]);
//        }
//        return out;
//    }
//
//    private static double[] getNegBin(double[] array) {
////        double[] returnArray = new double[2];
//        double mean = DiscreteStatistics.mean(array);
//        double variance = DiscreteStatistics.variance(array, mean);
//        double returnArray0 = (1 - (mean / variance));
//        double returnArray1 = (mean * ((1 - returnArray0) / returnArray0));
//        return new double[]{returnArray1, (returnArray0 / (1 - returnArray0))};
//    }
    public String output() {
        StringBuffer sb = new StringBuffer();
        sb.append("site\tcS\tuS\tcN\tuN\tsmooth(cS)\tsmooth(uS)\tsmooth(cN)\tsmooth(uN)" +
                "\tsmooth(cS/uS)\tsmooth(cN/uN)\tsmooth((cN/uN)/(cS/uS))\t[hpd]\n");
        for (int site = 0; site < numSites; site++) {
            sb.append(site + 1 + "\t" +
                      rawMeanStats[0][site] + "\t" + rawMeanStats[1][site] + "\t" + rawMeanStats[2][site] + "\t" + rawMeanStats[3][site] + "\t" +
                     smoothMeanStats[0][site] + "\t" + smoothMeanStats[1][site] + "\t" + smoothMeanStats[2][site] + "\t" + smoothMeanStats[3][site] + "\t" +
                     smoothMeanDnDsStats[0][site] + "\t" + smoothMeanDnDsStats[1][site] + "\t" + smoothMeanDnDsStats[2][site] + "\t" +
                     "["+smoothHPDDnDsStats[2][site][0] + "," + smoothHPDDnDsStats[2][site][1] + "]");
            if (smoothHPDDnDsStats[2][site][0] > 1 || smoothHPDDnDsStats[2][site][1] < 1) {
                sb.append("\t*");
            }
            sb.append("\n");
        }
        return sb.toString();
    }

    private static double[][] getArrayHPDintervals(double[][] array) {
        double[][] returnArray = new double[array.length][2];
        for (int row = 0; row < array.length; row++) {
            int counter = 0;
            for (int col = 0; col < array[0].length; col++) {
                if (!(((Double) array[row][col]).isNaN())) {
                    counter += 1;
                }
            }

            if (counter > 0) {
                double[] rowNoNaNArray = new double[counter];
                int index = 0;
                for (int col = 0; col < array[0].length; col++) {
                    if (!(((Double) array[row][col]).isNaN())) {
                        rowNoNaNArray[index] = array[row][col];
                        index += 1;
                    }
                }
                int[] indices = new int[counter];
                HeapSort.sort(rowNoNaNArray, indices);
                double hpdBinInterval[] = getHPDInterval(0.95, rowNoNaNArray, indices);

                returnArray[row][0] = hpdBinInterval[0];
                returnArray[row][1] = hpdBinInterval[1];
            } else {
                returnArray[row][0] = Double.NaN;
                returnArray[row][1] = Double.NaN;
            }

        }

        return returnArray;
    }

    private static double[] getHPDInterval(double proportion, double[] array, int[] indices) {

        double returnArray[] = new double[2];
        double minRange = Double.MAX_VALUE;
        int hpdIndex = 0;

        int diff = (int) Math.round(proportion * (double) array.length);
        for (int i = 0; i <= (array.length - diff); i++) {
            double minValue = array[indices[i]];
            double maxValue = array[indices[i + diff - 1]];
            double range = Math.abs(maxValue - minValue);
            if (range < minRange) {
                minRange = range;
                hpdIndex = i;
            }
        }
        returnArray[0] = array[indices[hpdIndex]];
        returnArray[1] = array[indices[hpdIndex + diff - 1]];
        return returnArray;
    }

    private static double[] mean(double[][] x) {
        double[] returnArray = new double[x.length];
        //System.out.println("lala");
        for (int i = 0; i < x.length; i++) {
            returnArray[i] = DiscreteStatistics.mean(x[i]);
            //System.out.println((i+1)+"\t"+DiscreteStatistics.mean(x[i]));
            //System.out.println(DiscreteStatistics.mean(x[i]));
        }
        return returnArray;
    }

    //   private static void print2DArray(double[][] array, String name) {
    //       try {
    //           PrintWriter outFile = new PrintWriter(new FileWriter(name), true);
    //
    //           for (int i = 0; i < array.length; i++) {
    //               for (int j = 0; j < array[0].length; j++) {
    //                   outFile.print(array[i][j]+"\t");
    //               }
    //          outFile.println("");
    //           }
    //           outFile.close();
    //
    //       } catch(IOException io) {
    //          System.err.print("Error writing to file: " + name);
    //       }
    //   }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return DNDS_PERSITE_ANALYSIS;
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
                long burnin = xo.getAttribute(MarginalLikelihoodAnalysisParser.BURN_IN, maxState / 10);
                //TODO: implement custom burn-in

                if (burnin < 0 || burnin >= maxState) {
                    burnin = maxState / 5;
                    System.out.println("WARNING: Burn-in larger than total number of states - using to 20%");
                }

                traces.setBurnIn(burnin);

                double samples[][][] = new double[NUM_VARIABLES][][];

                for (int variable = 0; variable < NUM_VARIABLES; ++variable) {
                    XMLObject cxo = xo.getChild(names[variable]);
                    String columnName = cxo.getStringAttribute(Attribute.NAME);
                    int traceStartIndex = -1;
                    int traceEndIndex = -1;
                    boolean traceIndexFound = false;

                    for (int i = 0; i < traces.getTraceCount(); i++) {
                        String traceName = traces.getTraceName(i);
                        if (traceName.trim().contains(columnName)) {
                            traceEndIndex = i;
                            if (!traceIndexFound) {
                                traceStartIndex = i;
                                traceIndexFound = true;
                            }
                        }
                    }
                    if (traceStartIndex == -1) {
                        throw new XMLParseException(columnName + " columns can not be found for " + getParserName() + " element.");
                    }
                    int numberOfSites = 1 + (traceEndIndex - traceStartIndex);

                    double[][] countPerSite = new double[numberOfSites][];
                    for (int a = 0; a < numberOfSites; a++) {
                        List<Double> values = (List)traces.getValues((a + traceStartIndex));
                        countPerSite[a] = new double[values.size()];
                        for (int i = 0; i < values.size(); i++) {
                            countPerSite[a][i] = values.get(i);
                        }
                    }

                    samples[variable] = countPerSite;

                }

                OldDnDsPerSiteAnalysis analysis = new OldDnDsPerSiteAnalysis(samples[COND_S], samples[UNCOND_S],
                        samples[COND_N], samples[UNCOND_N]);

                System.out.println(analysis.output());
                //TODO: save to file

                return analysis;

            } catch (FileNotFoundException fnfe) {
                throw new XMLParseException("File '" + fileName + "' can not be opened for " + getParserName() + " element.");
            } catch (IOException ioe) {
                throw new XMLParseException(ioe.getMessage());
            } catch (TraceException e) {
                throw new XMLParseException(e.getMessage());
            }
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "Performs a trace dN dS analysis.";
        }

        public Class getReturnType() {
            return OldDnDsPerSiteAnalysis.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                new StringAttributeRule(FileHelpers.FILE_NAME,
                        "The traceName of a BEAST log file (can not include trees, which should be logged separately"),
//                new ElementRule(UNCONDITIONAL_S_COLUMN, new XMLSyntaxRule[]{
//                       new StringAttributeRule(Attribute.NAME, "The column name")}),
//                new ElementRule(UNCONDITIONAL_N_COLUMN, new XMLSyntaxRule[]{
//                        new StringAttributeRule(Attribute.NAME, "The column name")}),
        };
    };

    final private int numSites;
    final private int numSamples;

    private double[][][] allSamples;
    final private static int NUM_VARIABLES = 4;
    final private static int COND_S = 0;
    final private static int UNCOND_S = 1;
    final private static int COND_N = 2;
    final private static int UNCOND_N = 3;
    final private static String[] names = {COND_SPERSITE_COLUMNS, UNCOND_SPERSITE_COLUMNS, COND_NPERSITE_COLUMNS, UNCOND_NPERSITE_COLUMNS};
    private double[][][] smoothSamples;
    private double[][][] smoothDnDsSamples;

    private static final boolean DEBUG = true;

    private double[][] rawMeanStats;
    private double[][] smoothMeanStats;
    private double[][] smoothMeanDnDsStats;
    private double[][][] smoothHPDDnDsStats;

}