/*
 * MarginalLikelihoodAnalysis.java
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

import dr.math.LogTricks;
import dr.math.MathUtils;
import dr.util.TaskListener;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Marc Suchard
 * @author Alexei Drummond
 *         <p/>
 *         Source translated from model_P.c (a component of BAli-Phy by Benjamin Redelings and Marc Suchard
 */
public class MarginalLikelihoodAnalysis {

    private final String traceName;
    private final List<Double> sample;
    private final long burnin;
    private final String analysisType; // "harmonic" for harmonic mean, "smoothed" for smoothed harmonic mean, "aicm" for AICM, "arithmetic" for arithmetic mean 
    private final int bootstrapLength;

    private boolean marginalLikelihoodCalculated = false;
    private double logMarginalLikelihood;
    private double bootstrappedSE;


//    public MarginalLikelihoodAnalysis(double[] sample, String traceName, int burnin) {
//        this(sample, traceName, burnin, false, 1000);
//    }

    public String getTraceName() {
        return traceName;
    }

    public long getBurnin() {
        return burnin;
    }

    /**
     * Constructor
     *
     * @param sample
     * @param traceName       used for 'toString' display purposes only
     * @param burnin          used for 'toString' display purposes only
     * @param analysisType
     * @param bootstrapLength a value of zero will turn off bootstrapping
     */
    public MarginalLikelihoodAnalysis(List<Double> sample, String traceName, long burnin, String analysisType, int bootstrapLength) {
        this.sample = sample;
        this.traceName = traceName;
        this.burnin = burnin;
        this.analysisType = analysisType;
        this.bootstrapLength = bootstrapLength;
//        System.err.println("setting burnin to "+burnin);
    }

    public double calculateLogMarginalLikelihood(List<Double> sample) {
        if (analysisType.equals("aicm")) {
            return logMarginalLikelihoodAICM(sample);
        } else if (analysisType.equals("smoothed")) {
            return logMarginalLikelihoodSmoothed(sample);
        } else if (analysisType.equals("arithmetic")) {
        	return logMarginalLikelihoodArithmetic(sample);
        } else {
            return logMarginalLikelihoodHarmonic(sample);
        }
    }
    
    /**
     * Calculates the log marginal likelihood of a model using the arithmetic mean estimator
     *
     * @param v a posterior sample of logLikelihoods
     * @return the log marginal likelihood
     */
    public double logMarginalLikelihoodArithmetic(List<Double> v) {
    	
    	int size = v.size();
    	
    	double sum = LogTricks.logZero;
        
        for (int i = 0; i < size; i++) {
            if (!Double.isNaN(v.get(i)) && !Double.isInfinite(v.get(i))) {
                sum = LogTricks.logSum(sum, v.get(i));
            } else {
                size--;
            }
        }
        
        return sum - StrictMath.log(size);
    	
    }

    /**
     * Calculates the log marginal likelihood of a model using Newton and Raftery's harmonic mean estimator
     *
     * @param v a posterior sample of logLikelihoods
     * @return the log marginal likelihood
     */

    public double logMarginalLikelihoodHarmonic(List<Double> v) {

        double sum = 0;
        final int size = v.size();
        for (int i = 0; i < size; i++)
            sum += v.get(i);

        double denominator = LogTricks.logZero;

        for (int i = 0; i < size; i++)
            denominator = LogTricks.logSum(denominator, sum - v.get(i));

        return sum - denominator + StrictMath.log(size);
    }

    /**
     * Calculates the AICM of a model using method-of-moments from Raftery et al. (2007)
     *
     * @param v a posterior sample of logLikelihoods
     * @return the AICM (lower values are better)
     */

    public double logMarginalLikelihoodAICM(List<Double> v) {

        double sum = 0;
        final int size = v.size();
        for (int i = 0; i < size; i++)
            sum += v.get(i);

        double mean = sum / (double) size;

        double var = 0;
        for (int i = 0; i < size; i++)
            var += (v.get(i)-mean) * (v.get(i)-mean);
        var /= (double) size - 1;

        return 2 * var - 2 * mean;

    }

    public void calculate() {

        logMarginalLikelihood = calculateLogMarginalLikelihood(sample);
        if (bootstrapLength > 1) {
            final int sampleLength = sample.size();
            List<Double> bsSample = new ArrayList<Double>();
            Double[] bootstrappedLogML = new Double[bootstrapLength];
            double sum = 0;

            double progress = 0.0;
            double delta = 1.0 / bootstrapLength;
            
            //System.err.println("HME = " + logMarginalLikelihood);

            for (int i = 0; i < bootstrapLength; i++) {
            //	if (i % 10 == 0) {
            //		System.err.println((i+1) + "/" + bootstrapLength);
            //	}
            	
                fireProgress(progress);
                progress += delta;

                int[] indices = MathUtils.sampleIndicesWithReplacement(sampleLength);
                for (int k = 0; k < sampleLength; k++) {
                    bsSample.add(sample.get(indices[k]));
                }
                bootstrappedLogML[i] = calculateLogMarginalLikelihood(bsSample);
                sum += bootstrappedLogML[i];
                bsSample.clear();
            }
            sum /= bootstrapLength;
            double bootstrappedAverage = sum;
            // Summarize bootstrappedLogML
            double var = 0;
            for (int i = 0; i < bootstrapLength; i++) {
                var += (bootstrappedLogML[i] - bootstrappedAverage) *
                        (bootstrappedLogML[i] - bootstrappedAverage);
            }
            var /= (bootstrapLength - 1.0);
            bootstrappedSE = Math.sqrt(var);
        }

        fireProgress(1.0);

        marginalLikelihoodCalculated = true;
    }

    /**
     * Calculates the log marginal likelihood of a model using Newton and Raftery's smoothed estimator
     *
     * @param v     a posterior sample of logLikelihood
     * @param delta proportion of pseudo-samples from the prior
     * @param Pdata current estimate of the log marginal likelihood
     * @return the log marginal likelihood
     */
    @SuppressWarnings({"SuspiciousNameCombination"})
    public double logMarginalLikelihoodSmoothed(List<Double> v, double delta, double Pdata) {

        final double logDelta = StrictMath.log(delta);
        final double logInvDelta = StrictMath.log(1.0 - delta);
        final int n = v.size();
        final double logN = StrictMath.log(n);

        final double offset = logInvDelta - Pdata;

        double bottom = logN + logDelta - logInvDelta;
        double top = bottom + Pdata;

        for (int i = 0; i < n; i++) {
            double weight = -LogTricks.logSum(logDelta, offset + v.get(i));
            top = LogTricks.logSum(top, weight + v.get(i));
            bottom = LogTricks.logSum(bottom, weight);
        }

        return top - bottom;
    }

    public double getLogMarginalLikelihood() {
        if (!marginalLikelihoodCalculated) {
            calculate();
        }
        return logMarginalLikelihood;
    }

    public double getBootstrappedSE() {
        if (!marginalLikelihoodCalculated) {
            calculate();
        }
        return bootstrappedSE;
    }

    public String toString() {
        StringBuilder sb = new StringBuilder();
        
        if (analysisType.equals("smoothed")) {
            sb.append("log marginal likelihood (using smoothed harmonic mean)");
        }
        else if (analysisType.equals("aicm")) {
            sb.append("AICM");
        }
        else if (analysisType.equals("arithmetic")) {
        	sb.append("log marginal likelihood (using arithmetic mean)");
        }
        else {
            sb.append("log marginal likelihood (using harmonic mean)");
        }
        
        sb.append(" from ")
                .append(traceName)
                .append(" = ")
                .append(String.format("%5.4f", getLogMarginalLikelihood()));
        if (bootstrapLength > 1) {
            sb.append(" +/- ")
                    .append(String.format("%5.4f", getBootstrappedSE()));
        } else {
            sb.append("           ");
        }

        sb.append(" burnin=").append(burnin);
        if (bootstrapLength > 1)
            sb.append(" replicates=").append(bootstrapLength);
//        sb.append("\n");

        return sb.toString();

    }

    public double logMarginalLikelihoodSmoothed(List<Double> v) {

        final double delta = 0.01;  // todo make class adjustable by accessor/setter

        // Start with harmonic estimator as first guess
        double Pdata = logMarginalLikelihoodHarmonic(v);

        double deltaP = 1.0;

        int iterations = 0;

        double dx;

        final double tolerance = 1E-3; // todo make class adjustable by accessor/setter

        while (Math.abs(deltaP) > tolerance) {
            double g1 = logMarginalLikelihoodSmoothed(v, delta, Pdata) - Pdata;
            double Pdata2 = Pdata + g1;
            dx = g1 * 10.0;
            double g2 = logMarginalLikelihoodSmoothed(v, delta, Pdata + dx) - (Pdata + dx);
            double dgdx = (g2 - g1) / dx; // find derivative at Pdata

            double Pdata3 = Pdata - g1 / dgdx; // find new evaluation point
            if (Pdata3 < 2.0 * Pdata || Pdata3 > 0 || Pdata3 > 0.5 * Pdata) // step is too large
                Pdata3 = Pdata + 10.0 * g1;

            double g3 = logMarginalLikelihoodSmoothed(v, delta, Pdata3) - Pdata3;

            // Try to do a Newton's method step
            if (Math.abs(g3) <= Math.abs(g2) && ((g3 > 0) || (Math.abs(dgdx) > 0.01))) {
                deltaP = Pdata3 - Pdata;
                Pdata = Pdata3;
            }  // otherwise try to go 10 times as far as one step
            else if (Math.abs(g2) <= Math.abs(g1)) {
                Pdata2 += g2;
                deltaP = Pdata2 - Pdata;
                Pdata = Pdata2;
            }  // otherwise go just one step
            else {
                deltaP = g1;
                Pdata += g1;
            }

            iterations++;

            if (iterations > 400) { // todo make class adjustable by acessor/setter
                System.err.println("Probabilities are not converging!!!"); // todo should throw exception
                return LogTricks.logZero;
            }
        }
        return Pdata;
    }

    private TaskListener listener = null;

    public void setTaskListener(TaskListener listener) {
        this.listener = listener;
    }

    private void fireProgress(double progress) {
        if (listener != null) {
            listener.progress(progress);
        }
    }

}
