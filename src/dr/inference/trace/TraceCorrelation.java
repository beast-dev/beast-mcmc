/*
 * TraceCorrelation.java
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

import java.util.List;
import java.util.Map;

/**
 * A class that stores the correlation statistics for a trace.
 * The difference to TraceDistribution is mainly to add ACT and ESS
 * which require stepSize.
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: TraceCorrelation.java,v 1.2 2006/11/29 14:53:53 rambaut Exp $
 */
public class TraceCorrelation extends TraceDistribution {
    private final long stepSize;

    public TraceCorrelation(List<Double> values, TraceType traceType, long stepSize) {
        super(values, traceType);
        this.stepSize = stepSize;

        if (stepSize > 0) {
            analyseCorrelation(values, stepSize);
        }
    }

    public TraceCorrelation(List<Double> values, Map<Integer, String> categoryLabelMap, List<Integer> categoryOrder, long stepSize) {
        super(values, categoryLabelMap, categoryOrder);

        this.stepSize = stepSize;

        if (stepSize > 0) {
            analyseCorrelation(values, stepSize);
        }
    }

    public double getStdErrorOfMean() {
        return stdErrorOfMean;
    }

    public double getACT() {
        return ACT;
    }

    public double getESS() {
        return ESS;
    }

    //************************************************************************
    // private methods
    //************************************************************************

    private double stdErrorOfMean;
    private double stdErrorOfVariance;
    private double ACT;
    private double stdErrOfACT;
    private double ESS;

    private static final int MAX_LAG = 2000;

    private void analyseCorrelation(List<Double> values, long stepSize) {
//        this.values = values; // move to TraceDistribution(T[] values)

        if (stepSize > 0) {
            if (getTraceType().isNumber()) {
                double[] doubleValues = new double[values.size()];
                for (int i = 0; i < values.size(); i++) {
                    doubleValues[i] = values.get(i);
                }
                analyseCorrelationNumeric(doubleValues, stepSize);

            } else if (getTraceType() == TraceType.CATEGORICAL) {
                //todo Do not know how to calculate
                stdErrorOfMean = Double.NaN;
                ACT = Double.NaN;
                ESS = Double.NaN;
                stdErrOfACT = Double.NaN;
//            throw new UnsupportedOperationException("should not be categorical");
            } else {
                throw new RuntimeException("Trace type is not recognized");
            }
        }
    }

    /**
     * Analyze trace for numeric values
     *
     * @param values   the values
     * @param stepSize the sampling frequency of the values
     */
    private void analyseCorrelationNumeric(double[] values, long stepSize) {

        final int samples = values.length;
        int maxLag = Math.min(samples - 1, MAX_LAG);

        double[] gammaStat = new double[maxLag];
        //double[] varGammaStat = new double[maxLag];
        double varStat = 0.0;
        //double varVarStat = 0.0;
        //double assVarCor = 1.0;
        //double del1, del2;

        for (int lag = 0; lag < maxLag; lag++) {
            for (int j = 0; j < samples - lag; j++) {
                final double del1 = values[j] - mean;
                final double del2 = values[j + lag] - mean;
                gammaStat[lag] += (del1 * del2);
                //varGammaStat[lag] += (del1*del1*del2*del2);
            }

            gammaStat[lag] /= ((double) (samples - lag));
            //varGammaStat[lag] /= ((double) samples-lag);
            //varGammaStat[lag] -= (gammaStat[0] * gammaStat[0]);

            if (lag == 0) {
                varStat = gammaStat[0];
                //varVarStat = varGammaStat[0];
                //assVarCor = 1.0;
            } else if (lag % 2 == 0) {
                // fancy stopping criterion :)
                if (gammaStat[lag - 1] + gammaStat[lag] > 0) {
                    varStat += 2.0 * (gammaStat[lag - 1] + gammaStat[lag]);
                    // varVarStat += 2.0*(varGammaStat[lag-1] + varGammaStat[lag]);
                    // assVarCor  += 2.0*((gammaStat[lag-1] * gammaStat[lag-1]) + (gammaStat[lag] * gammaStat[lag])) / (gammaStat[0] * gammaStat[0]);
                }
                // stop
                else {
                    maxLag = lag;
                }
            }
        }

        // standard error of mean
        stdErrorOfMean = Math.sqrt(varStat / samples);

        // auto correlation time
        if (gammaStat[0] == 0) {
            ACT = 0;
        } else {
            ACT = stepSize * varStat / gammaStat[0];
        }

        // effective sample size
        if (ACT == 0) {
            ESS = 1;
        } else {
            ESS = (stepSize * samples) / ACT;
        }

        // standard deviation of autocorrelation time
        stdErrOfACT = (2.0 * Math.sqrt(2.0 * (2.0 * (double) (maxLag + 1)) / samples) * (varStat / gammaStat[0]) * stepSize);

//        minEqualToMax = true;
    }

}