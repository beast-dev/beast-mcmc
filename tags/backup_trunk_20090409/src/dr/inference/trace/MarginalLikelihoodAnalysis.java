/*
 * MarginalLikelihoodAnalysis.java
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

import dr.math.LogTricks;
import dr.math.MathUtils;
import dr.util.Attribute;
import dr.util.TaskListener;
import dr.xml.*;

import java.io.File;
import java.io.FileNotFoundException;

/**
 * @author Marc Suchard
 * @author Alexei Drummond
 *         <p/>
 *         Source translated from model_P.c (a component of BAli-Phy by Benjamin Redelings and Marc Suchard
 */
public class MarginalLikelihoodAnalysis {

    public static final String ML_ANALYSIS = "marginalLikelihoodAnalysis";
    public static final String FILE_NAME = "fileName";
    public static final String BURN_IN = "burnIn";
    public static final String COLUMN_NAME = "likelihoodColumn";
    public static final String DO_BOOTSTRAP = "bootstrap";
    public static final String ONLY_HARMONIC = "harmonicOnly";
    public static final String BOOTSTRAP_LENGTH = "bootstrapLength";

    private String traceName;
    private double[] sample;
    private int burnin;
    private boolean harmonicOnly;
    private int bootstrapLength;


    private boolean marginalLikelihoodCalculated = false;
    private double logMarginalLikelihood;
    private double bootstrappedSE;


    public MarginalLikelihoodAnalysis(double[] sample, String traceName, int burnin) {
        this(sample, traceName, burnin, false, 1000);
    }

    public String getTraceName() {
        return traceName;
    }

    public int getBurnin() {
        return burnin;
    }

    /**
     * Constructor
     *
     * @param sample
     * @param traceName       used for 'toString' display purposes only
     * @param burnin          used for 'toString' display purposes only
     * @param harmonicOnly
     * @param bootstrapLength a value of zero will turn off bootstrapping
     */
    public MarginalLikelihoodAnalysis(double[] sample, String traceName, int burnin, boolean harmonicOnly, int bootstrapLength) {
        this.sample = sample;
        this.traceName = traceName;
        this.burnin = burnin;
        this.harmonicOnly = harmonicOnly;
        this.bootstrapLength = bootstrapLength;
//        System.err.println("setting burnin to "+burnin);
    }

    public double calculateLogMarginalLikelihood(double[] sample) {
        if (harmonicOnly)
            return logMarginalLikelihoodHarmonic(sample);
        else
            return logMarginalLikelihoodSmoothed(sample);
    }

    /**
     * Calculates the log marginal likelihood of a model using Newton and Raftery's harmonic mean estimator
     *
     * @param v a posterior sample of logLikelihoods
     * @return the log marginal likelihood
     */

    public double logMarginalLikelihoodHarmonic(double[] v) {

        double sum = 0;
        final int size = v.length;
        for (int i = 0; i < size; i++)
            sum += v[i];

        double denominator = LogTricks.logZero;

        for (int i = 0; i < size; i++)
            denominator = LogTricks.logSum(denominator, sum - v[i]);

        return sum - denominator + StrictMath.log(size);
    }

    public void calculate() {

        logMarginalLikelihood = calculateLogMarginalLikelihood(sample);
        if (bootstrapLength > 1) {
            final int sampleLength = sample.length;
            double[] bsSample = new double[sampleLength];
            double[] bootstrappedLogML = new double[bootstrapLength];
            double sum = 0;

            double progress = 0.0;
            double delta = 1.0 / bootstrapLength;

            for (int i = 0; i < bootstrapLength; i++) {
                fireProgress(progress);
                progress += delta;

                int[] indices = MathUtils.sampleIndicesWithReplacement(sampleLength);
                for (int k = 0; k < sampleLength; k++)
                    bsSample[k] = sample[indices[k]];
                bootstrappedLogML[i] = calculateLogMarginalLikelihood(bsSample);
                sum += bootstrappedLogML[i];
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
     * @param v     a posterior sample of logLilelihood
     * @param delta proportion of pseudo-samples from the prior
     * @param Pdata current estimate of the log marginal likelihood
     * @return the log marginal likelihood
     */
    @SuppressWarnings({"SuspiciousNameCombination"})
    public double logMarginalLikelihoodSmoothed(double[] v, double delta, double Pdata) {

        final double logDelta = StrictMath.log(delta);
        final double logInvDelta = StrictMath.log(1.0 - delta);
        final int n = v.length;
        final double logN = StrictMath.log(n);

        final double offset = logInvDelta - Pdata;

        double bottom = logN + logDelta - logInvDelta;
        double top = bottom + Pdata;

        for (int i = 0; i < n; i++) {
            double weight = -LogTricks.logSum(logDelta, offset + v[i]);
            top = LogTricks.logSum(top, weight + v[i]);
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
        sb.append("log P(")
                .append(traceName)
                .append("|Data) = ")
                .append(String.format("%5.4f", getLogMarginalLikelihood()));
        if (bootstrapLength > 1) {
            sb.append(" +/- ")
                    .append(String.format("%5.4f", getBootstrappedSE()));
        } else {
            sb.append("           ");
        }
        if (harmonicOnly)
            sb.append(" (harmonic)");
        else
            sb.append(" (smoothed)");
        sb.append(" burnin=").append(burnin);
        if (bootstrapLength > 1)
            sb.append(" replicates=").append(bootstrapLength);
//        sb.append("\n");

        return sb.toString();

    }

    public double logMarginalLikelihoodSmoothed(double[] v) {

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

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return ML_ANALYSIS;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            String fileName = xo.getStringAttribute(FILE_NAME);
            try {

                File file = new File(fileName);
                String name = file.getName();
                String parent = file.getParent();

                if (!file.isAbsolute()) {
                    parent = System.getProperty("user.dir");
                }

                file = new File(parent, name);

                fileName = file.getAbsolutePath();

                XMLObject cxo = (XMLObject) xo.getChild(COLUMN_NAME);
                String likelihoodName = cxo.getStringAttribute(Attribute.NAME);

                LogFileTraces traces = new LogFileTraces(fileName, file);
                traces.loadTraces();
                int maxState = traces.getMaxState();

                // leaving the burnin attribute off will result in 10% being used
                int burnin = xo.getAttribute(BURN_IN, maxState / 10);

                if (burnin < 0 || burnin >= maxState) {
                    burnin = maxState / 10;
                    System.out.println("WARNING: Burn-in larger than total number of states - using to 10%");
                }

                traces.setBurnIn(burnin);

                int traceIndex = -1;
                for (int i = 0; i < traces.getTraceCount(); i++) {
                    String traceName = traces.getTraceName(i);
                    if (traceName.equals(likelihoodName)) {
                        traceIndex = i;
                        break;
                    }
                }

                if (traceIndex == -1) {
                    throw new XMLParseException("Column '" + likelihoodName + "' can not be found for " + getParserName() + " element.");
                }

                boolean harmonicOnly = false;
                if (cxo.hasAttribute(ONLY_HARMONIC)) {
                    harmonicOnly = cxo.getBooleanAttribute(ONLY_HARMONIC);
                }

                int bootstrapLength = cxo.getAttribute(BOOTSTRAP_LENGTH, 1000);

                double sample[] = new double[traces.getStateCount()];
                traces.getValues(traceIndex, sample);

                MarginalLikelihoodAnalysis analysis = new MarginalLikelihoodAnalysis(
                        sample,
                        traces.getTraceName(traceIndex), burnin,
                        harmonicOnly, bootstrapLength);

                System.out.println(analysis.toString());

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
            return "Performs a trace analysis. Estimates the mean of the various statistics in the given log file.";
        }

        public Class getReturnType() {
            return MarginalLikelihoodAnalysis.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                new StringAttributeRule(FILE_NAME, "The traceName of a BEAST log file (can not include trees, which should be logged separately"),
                AttributeRule.newIntegerRule("burnIn", true)
                //, "The number of states (not sampled states, but actual states) that are discarded from the beginning of the trace before doing the analysis" ),
		};
	};
}
