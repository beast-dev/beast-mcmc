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

import dr.stats.DiscreteStatistics;
import dr.util.HeapSort;
import dr.util.NumberFormatter;

import java.io.FileReader;
import java.io.Reader;

/**
 * @author Alexei Drummond
 *
 * @version $Id: TraceAnalysis.java,v 1.23 2005/05/24 20:26:00 rambaut Exp $
 */
public class TraceAnalysis {

    public TraceAnalysis(Trace trace, int burnIn) {
        this(trace, burnIn, true);
    }

    public TraceAnalysis(Trace trace, int burnin, boolean correlation) {
        this.correlation = correlation;
        this.trace = trace;
        statName = trace.getId();
        maxState = Trace.Utils.getMaximumState(trace);
        stepSize = trace.getStepSize();

        if (burnin < 0 || burnin >= maxState) {
            burnin = maxState / 10;
            System.out.println("WARNING: Burn-in larger than total number of states - using to 10%");
        }

        reanalyze(burnin);
    }

    /**
     * Actually analyzes the trace given the burnin
     */
    public void reanalyze(int burnin) {

        if (burnin >= 0 && burnin < maxState) {
            this.burnin = burnin;

            double[] values = trace.getValues(burnin);
            mean = DiscreteStatistics.mean(values);

            int[] indices = new int[values.length];
            HeapSort.sort(values, indices);
            median = DiscreteStatistics.quantile(0.5, values, indices);
            cpdLower = DiscreteStatistics.quantile(0.025, values, indices);
            cpdUpper = DiscreteStatistics.quantile(0.975, values, indices);
            calculateHPDInterval(0.95, values, indices);

            if (correlation) analyze(values, trace.getStepSize());
        }
    }

    public int getBurnin() { return burnin; }
    public double getHPDUpper() { return hpdUpper; }
    public double getHPDLower() { return hpdLower; }
    public double getCPDUpper() { return cpdUpper; }
    public double getCPDLower() { return cpdLower; }
    public double getMean() { return mean; }
    public double getMedian() { return median; }

    public double getStdError() {
        if (!correlation) throw new RuntimeException("Correlation-based statistics not calculated!");
        return stdErrorOfMean;
    }

    public double getAutoCorrelationTime() {
        if (!correlation) throw new RuntimeException("Correlation-based statistics not calculated!");
        return ACT;
    }
    public double getEffectiveSampleSize() {
        if (!correlation) throw new RuntimeException("Correlation-based statistics not calculated!");
        return ESS;
    }

    public String getName() { return statName; }
    public int getMaximumState() { return maxState; }
    public int getStepSize() { return stepSize; }

    public Trace getTrace() { return trace; }

    /**
     * @return an array og analyses of the statistics in a log file.
     */
    public static TraceAnalysis[] analyzeLogFile(Reader reader, int burnin) throws java.io.IOException {
        return analyzeLogFile(reader, burnin, false);
    }

    /**
     * @return an array og analyses of the statistics in a log file.
     */
    public static TraceAnalysis[] analyzeLogFile(Reader reader, int burnin, boolean percentage) throws java.io.IOException {

        Trace[] traces = Trace.Utils.loadTraces(reader);
        reader.close();

        int maxState = Trace.Utils.getMaximumState(traces[0]);

        if (burnin == -1) {
            burnin = maxState / 10;
        } else if (percentage) {
            burnin = (int)Math.round((double)maxState * (double)burnin / 100.0);
        }

        if (burnin < 0 || burnin >= maxState) {
            burnin = maxState / 10;
            System.out.println("WARNING: Burn-in larger than total number of states - using to 10%");
        }

        TraceAnalysis[] analysis = new TraceAnalysis[traces.length];

        for (int j = 0; j < traces.length; j++) {
            analysis[j] = new TraceAnalysis(traces[j], burnin);
        }

        return analysis;
    }

    /**
     * @return an array og analyses of the statistics in a log file.
     */
    public static TraceAnalysis[] analyzeLogFile(String fileName, int burnin) throws java.io.IOException {

        FileReader reader = new FileReader(fileName);
        return analyzeLogFile(reader, burnin);
    }

    public static TraceAnalysis[] report(String fileName) throws java.io.IOException {

        return report(new FileReader(fileName), -1);
    }

    public static TraceAnalysis[] report(String fileName, int burnin) throws java.io.IOException {

        return report(new FileReader(fileName), burnin);
    }

    /**
     * @param reader the log file to report
     * @param burnin the number of states of burnin or if -1 then use 10%
     */
    public static TraceAnalysis[] report(Reader reader, int burnin) throws java.io.IOException {

        int fieldWidth = 14;
        int firstField = 25;
        NumberFormatter formatter = new NumberFormatter(6);
        formatter.setPadding(true);
        formatter.setFieldWidth(fieldWidth);

        TraceAnalysis[] analysis = analyzeLogFile(reader, burnin);

        int maxState = analysis[0].maxState;
        burnin = analysis[0].burnin;

        System.out.println();
        System.out.println("burnIn=" + burnin);
        System.out.println("maxState=" + maxState);
        System.out.println();

        System.out.print(formatter.formatToFieldWidth("statistic", firstField));
        String[] names = new String[] {"mean", "hpdLower", "hpdUpper", "ESS"};

        for (int i =0; i < names.length; i++) {
            System.out.print(formatter.formatToFieldWidth(names[i], fieldWidth));
        }
        System.out.println();

        int warning = 0;
        for (int j = 0; j < analysis.length; j++) {
            System.out.print(formatter.formatToFieldWidth(analysis[j].getName(), firstField));
            System.out.print(formatter.format(analysis[j].getMean()));
            System.out.print(formatter.format(analysis[j].getHPDLower()));
            System.out.print(formatter.format(analysis[j].getHPDUpper()));
            System.out.print(formatter.format(analysis[j].getEffectiveSampleSize()));
            double ess = analysis[j].getEffectiveSampleSize();
            if (ess < 100) {
                warning += 1;
                System.out.println("*");
            } else {
                System.out.println();
            }
        }
        System.out.println();

        /*
          System.out.println("Correlation matrix");
          System.out.print(formatter.formatToFieldWidth("", firstField));
          for (int i =0; i < analysis.length; i++) {
              System.out.print(formatter.formatToFieldWidth(analysis[i].getName(), fieldWidth));
          }
          System.out.println();
          for (int i = 0; i < analysis.length; i++) {

              double[] x = analysis[i].getTrace().getValues(analysis[i].getBurnin());


              System.out.print(formatter.formatToFieldWidth(analysis[i].getName(), firstField));
              for (int j = 0; j < i; j++) {

                  double[] y = analysis[j].getTrace().getValues(analysis[j].getBurnin());
                  System.out.print(formatter.format(DiscreteStatistics.covariance(x, y)));
              }
              System.out.println();
          }
          */

        System.out.println();

        if (warning > 0) {
            System.out.println(" * WARNING: The results of this MCMC analysis may be invalid as ");
            System.out.println("            one or more statistics had very low effective sample sizes (ESS)");
        }

        System.out.flush();
        return analysis;
    }

    /**
     * @param reader the log file to report
     * @param burnin the number of states of burnin or if -1 then use 10%
     */
    public static TraceAnalysis[] shortReport(String name, Reader reader, int burnin, boolean drawHeader) throws java.io.IOException {

        TraceAnalysis[] analysis = analyzeLogFile(reader, burnin);

        int maxState = analysis[0].maxState;
        burnin = analysis[0].burnin;

        double minESS = Double.MAX_VALUE;

        if (drawHeader) {
            System.out.print("file\t");
            for (int j = 0; j < analysis.length; j++) {
                System.out.print(analysis[j].getName()+"\t");
            }
            System.out.println("minESS\tchainLength");
        }

        System.out.print(name + "\t");
        for (int j = 0; j < analysis.length; j++) {
            System.out.print(analysis[j].getMean()+"\t");
            double ess = analysis[j].getEffectiveSampleSize();
            if (ess < minESS) {
                minESS = ess;
            }
        }
        System.out.println(minESS+ "\t" + maxState);
        return analysis;
    }

    /**
     * @param proportion the proportion of probability mass oncluded within interval.
     */
    private void calculateHPDInterval(double proportion, double[] array, int[] indices) {

        double minRange = Double.MAX_VALUE;
        int hpdIndex = 0;

        int diff = (int)Math.round(proportion * (double)array.length);
        for (int i =0; i <= (array.length - diff); i++) {
            double minValue = array[indices[i]];
            double maxValue = array[indices[i+diff-1]];
            double range = Math.abs(maxValue - minValue);
            if (range < minRange) {
                minRange = range;
                hpdIndex = i;
            }
        }
        hpdLower = array[indices[hpdIndex]];
        hpdUpper = array[indices[hpdIndex+diff-1]];
    }

    /**
     * Analyze trace
     * @param statistic an array of the statistic
     * @param update the step size between instances of statistic
     */
    private void analyze(double[] statistic, int update) {

        int maxLag = MAX_OFFSET;
        int samples = statistic.length;
        if (samples-1 < maxLag) {
            maxLag = samples-1;
        }

        double[] gammaStat = new double[maxLag];
        double[] varGammaStat = new double[maxLag];
        double meanStat = 0.0;
        double varStat = 0.0;
        double varVarStat = 0.0;
        double assVarCor = 0.0;
        double del1, del2;

        for (int i = 0; i < samples; i++) {
            meanStat += statistic[i];
        }
        meanStat /= samples;

        for (int lag=0; lag < maxLag; lag++) {
            for (int j = 0; j < samples-lag; j++) {
                del1=statistic[j] - meanStat;
                del2=statistic[j+lag] - meanStat;
                gammaStat[lag] += ( del1*del2 );
                varGammaStat[lag] += (del1*del1*del2*del2);
            }

            gammaStat[lag] /= ((double)samples);
            varGammaStat[lag] /= ((double) samples-lag);
            varGammaStat[lag] -= (gammaStat[0] * gammaStat[0]);
	        
            if (lag==0) {
                varStat = gammaStat[0];
                varVarStat = varGammaStat[0];
                assVarCor = 1.0;
            }
            else if (lag%2==0)
            {
                // fancy stopping criterion :)
                if (gammaStat[lag-1] + gammaStat[lag] > 0) {
                    varStat    += 2.0*(gammaStat[lag-1] + gammaStat[lag]);
                    varVarStat += 2.0*(varGammaStat[lag-1] + varGammaStat[lag]);
                    assVarCor  += 2.0*((gammaStat[lag-1] * gammaStat[lag-1]) + (gammaStat[lag] * gammaStat[lag])) / (gammaStat[0] * gammaStat[0]);
                }
                // stop
                else
                    maxLag=lag;
            }
        }

        // standard error of mean
        stdErrorOfMean = Math.sqrt(varStat/samples);

        // variance of statistic
        //variance = gammaStat[0];

        // standard error of variance
        //stdErrorOfVariance = Math.sqrt(varVarStat/samples);

        // auto correlation time
        ACT = update * varStat / gammaStat[0];

        // effective sample size
        ESS = (update * statistic.length) / ACT;

        // M
        //M = lag;

        // M(updates)
        //M_update = lag * update;

        // standard deviation of autocorrelation time
        //stdErrOfACT = (2.0* Math.sqrt(2.0*(2.0*(double) maxLag+1)/samples)*(varStat/gammaStat[0])*update);

        //assymptotic std of correlation function
        //assStdOfCorrelationFunction = Math.sqrt(2.0*assVarCor/samples);
    }

    private double stdErrorOfMean;
    private double mean, median;
    //	private double variance;
    private int burnin = -1;
    private double ACT, ESS;
//    private double stdErrorOfVariance;
//	private double M, M_update;
//	private double stdErrOfACT;
    //	private double assStdOfCorrelationFunction;
    private double cpdLower, cpdUpper, hpdLower, hpdUpper;
    private static final int MAX_OFFSET = 2000;
    private boolean correlation;
    private String statName = null;
    private int maxState;
    private int stepSize;
    private Trace trace;
}