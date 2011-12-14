/*
 * TraceAnalysis.java
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

import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.List;

/**
 * @author Alexei Drummond
 * @version $Id: TraceAnalysis.java,v 1.23 2005/05/24 20:26:00 rambaut Exp $
 */
public class TraceAnalysis {

    /**
     * @param fileName the name of the log file to analyze
     * @param burnin   the state to discard up to
     * @return an array og analyses of the statistics in a log file.
     * @throws java.io.IOException if general error reading file
     * @throws TraceException      if trace file in wrong format or corrupted
     */
    public static LogFileTraces analyzeLogFile(String fileName, int burnin) throws java.io.IOException, TraceException {

        File file = new File(fileName);
        LogFileTraces traces = new LogFileTraces(fileName, file);
        traces.loadTraces();
        traces.setBurnIn(burnin);

        for (int i = 0; i < traces.getTraceCount(); i++) {
            traces.analyseTrace(i);
        }
        return traces;
    }

    public static TraceList report(String fileName) throws java.io.IOException, TraceException {
        return report(fileName, -1, null);
    }

    public static TraceList report(String fileName, int burnin, String likelihoodName) throws java.io.IOException, TraceException {
        return report(fileName, burnin, likelihoodName, true);
    }

    public static TraceList report(String fileName, int inBurnin, String likelihoodName, boolean withStdError)
            throws java.io.IOException, TraceException {

//        int fieldWidth = 14;
//        int firstField = 25;
//        NumberFormatter formatter = new NumberFormatter(4);
//        formatter.setPadding(true);
//        formatter.setFieldWidth(fieldWidth);

        File file = new File(fileName);

        LogFileTraces traces = new LogFileTraces(fileName, file);
//        if (traces == null) {
//            throw new TraceException("Trace file is empty.");
//        }
        traces.loadTraces();

//        traces.addTrace("R0", traces.getTraceIndex("bdss.psi"));

        int burnin = inBurnin;
        if (burnin == -1) {
            burnin = traces.getMaxState() / 10;
        }

        traces.setBurnIn(burnin);

//        System.out.println();
        System.out.println("burnIn   <= " + burnin + ",   maxState  = " + traces.getMaxState());
//        System.out.println();

        System.out.print("statistic");
        String[] names;

        if (!withStdError)
            names = new String[]{"mean", "hpdLower", "hpdUpper", "ESS"};
        else
            names = new String[]{"mean", "stdErr", "median", "hpdLower", "hpdUpper", "ESS", "50hpdLower", "50hpdUpper"};

        for (String name : names) {
            System.out.print("\t" + name);
        }
        System.out.println();

        int warning = 0;
        for (int i = 0; i < traces.getTraceCount(); i++) {
            traces.analyseTrace(i);
            TraceDistribution distribution = traces.getDistributionStatistics(i);

            double ess = distribution.getESS();
            System.out.print(traces.getTraceName(i));
            System.out.print("\t" + formattedNumber(distribution.getMean()));

            if (withStdError) {
                System.out.print("\t" + formattedNumber(distribution.getStdError()));
                System.out.print("\t" + formattedNumber(distribution.getMedian()));
            }

            System.out.print("\t" + formattedNumber(distribution.getLowerHPD()));
            System.out.print("\t" + formattedNumber(distribution.getUpperHPD()));

            System.out.print("\t" + formattedNumber(ess));
            
            if (withStdError) {
                System.out.print("\t" + formattedNumber(distribution.getHpdLowerCustom()));
                System.out.print("\t" + formattedNumber(distribution.getHpdUpperCustom()));
            }

            if (ess < 100) {
                warning += 1;
                System.out.println("\t" + "*");
            } else {
                System.out.println("\t");
            }
        }
        System.out.println();

        if (warning > 0) {
            System.out.println(" * WARNING: The results of this MCMC analysis may be invalid as ");
            System.out.println("            one or more statistics had very low effective sample sizes (ESS)");
        }

        if (likelihoodName != null) {
            System.out.println();
            int traceIndex = -1;
            for (int i = 0; i < traces.getTraceCount(); i++) {
                String traceName = traces.getTraceName(i);
                if (traceName.equals(likelihoodName)) {
                    traceIndex = i;
                    break;
                }
            }

            if (traceIndex == -1) {
                throw new TraceException("Column '" + likelihoodName +
                        "' can not be found for marginal likelihood analysis.");
            }

            String analysisType = "aicm";
            int bootstrapLength = 1000;

            List<Double> sample = traces.getValues(traceIndex);

            MarginalLikelihoodAnalysis analysis = new MarginalLikelihoodAnalysis(sample,
                    traces.getTraceName(traceIndex), burnin, analysisType, bootstrapLength);

            System.out.println(analysis.toString());
        }

        System.out.flush();
        return traces;
    }

    public static void reportTrace(String fileName, int inBurnin, String traceName) throws IOException, TraceException {
        File file = new File(fileName);

        LogFileTraces traces = new LogFileTraces(fileName, file);
        traces.loadTraces();
        int burnin = inBurnin;
        if (burnin == -1) {
            burnin = traces.getMaxState() / 10;
        }

        traces.setBurnIn(burnin);

//        System.out.println();
//        System.out.println("burnIn   <= " + burnin + ",   maxState  = " + traces.getMaxState());
//        System.out.println();

//        System.out.print("statistic");
//        String[] names = new String[]{"mean", "stdErr", "median", "hpdLower", "hpdUpper", "50hpdLower", "50hpdUpper"};//, "ESS"};
//
//        for (String name : names) {
//            System.out.print("\t" + name);
//        }
//        System.out.println();

        int id = traces.getTraceIndex(traceName);

        traces.analyseTrace(id);
        TraceDistribution distribution = traces.getDistributionStatistics(id);

        double ess = distribution.getESS();
//            System.out.print(traces.getTraceName(id) + "\t");
        System.out.print(formattedNumber(distribution.getMean()) + "\t");


        System.out.print(formattedNumber(distribution.getStdError()) + "\t");
        System.out.print(formattedNumber(distribution.getMedian()) + "\t");


        System.out.print(formattedNumber(distribution.getLowerHPD()) + "\t");
        System.out.print(formattedNumber(distribution.getUpperHPD()) + "\t");

        System.out.print(formattedNumber(distribution.getHpdLowerCustom()) + "\t");
        System.out.print(formattedNumber(distribution.getHpdUpperCustom()) + "\t");
        System.out.println();

//            System.out.print(SummaryStatisticsPanel.formattedNumber(ess));

    }

    /**
     * @param burnin         the number of states of burnin or if -1 then use 10%
     * @param filename       the file name of the log file to report on
     * @param drawHeader     if true then draw header
     * @param stdErr         if true then report the standard deviation of the mean
     * @param hpds           if true then report 95% hpd upper and lower
     * @param individualESSs minimum number of ESS with which to throw warning
     * @param likelihoodName column name
     * @return the traces loaded from given file to create this short report
     * @throws java.io.IOException if general error reading file
     * @throws TraceException      if trace file in wrong format or corrupted
     */
    public static TraceList shortReport(String filename,
                                        final int burnin, boolean drawHeader,
                                        boolean hpds, boolean individualESSs, boolean stdErr,
                                        String likelihoodName) throws java.io.IOException, TraceException {

        TraceList traces = analyzeLogFile(filename, burnin);

        int maxState = traces.getMaxState();

        double minESS = Double.MAX_VALUE;

        if (drawHeader) {
            System.out.print("file\t");
            for (int i = 0; i < traces.getTraceCount(); i++) {
                String traceName = traces.getTraceName(i);
                System.out.print(traceName + "\t");
                if (stdErr)
                    System.out.print(traceName + " stdErr\t");
                if (hpds) {
                    System.out.print(traceName + " hpdLower\t");
                    System.out.print(traceName + " hpdUpper\t");
                }
                if (individualESSs) {
                    System.out.print(traceName + " ESS\t");
                }
            }
            System.out.print("minESS\t");
            if (likelihoodName != null) {
                System.out.print("marginal likelihood\t");
                System.out.print("stdErr\t");
            }
            System.out.println("chainLength");
        }

        System.out.print(filename + "\t");
        for (int i = 0; i < traces.getTraceCount(); i++) {
            //TraceDistribution distribution = traces.getDistributionStatistics(i);
            TraceCorrelation distribution = traces.getCorrelationStatistics(i);
            System.out.print(distribution.getMean() + "\t");
            if (stdErr)
                System.out.print(distribution.getStdErrorOfMean() + "\t");
            if (hpds) {
                System.out.print(distribution.getLowerHPD() + "\t");
                System.out.print(distribution.getUpperHPD() + "\t");
            }
            if (individualESSs) {
                System.out.print(distribution.getESS() + "\t");
            }
            double ess = distribution.getESS();
            if (ess < minESS) {
                minESS = ess;
            }
        }

        System.out.print(minESS + "\t");

        if (likelihoodName != null) {
            int traceIndex = -1;
            for (int i = 0; i < traces.getTraceCount(); i++) {
                String traceName = traces.getTraceName(i);
                if (traceName.equals(likelihoodName)) {
                    traceIndex = i;
                    break;
                }
            }

            if (traceIndex == -1) {
                throw new TraceException("Column '" + likelihoodName + "' can not be found in file " + filename + ".");
            }

            String analysisType = "aicm";
            int bootstrapLength = 1000;

            List<Double> sample = traces.getValues(traceIndex);

            MarginalLikelihoodAnalysis analysis = new MarginalLikelihoodAnalysis(sample,
                    traces.getTraceName(traceIndex), burnin, analysisType, bootstrapLength);

            System.out.print(analysis.getLogMarginalLikelihood() + "\t");
            System.out.print(analysis.getBootstrappedSE() + "\t");
        }

        System.out.println(maxState);
        return traces;
    }

    public static String formattedNumber(double value) {
        DecimalFormat formatter = new DecimalFormat("0.####E0");
        DecimalFormat formatter2 = new DecimalFormat("####0.####");

        if (value > 0 && (Math.abs(value) < 0.01 || Math.abs(value) >= 100000.0)) {
            return formatter.format(value);
        } else return formatter2.format(value);
    }
}
