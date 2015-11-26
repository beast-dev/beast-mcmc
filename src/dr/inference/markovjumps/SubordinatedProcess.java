/*
 * SubordinatedProcess.java
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

package dr.inference.markovjumps;

import dr.math.GammaFunction;
import dr.math.MathUtils;
import dr.math.matrixAlgebra.Vector;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;


/**
 * A class to represent a Poisson process and discrete-time Markov chain that subordinate
 * a continuous-time Markov chain in the interval [0,T]. The subordinator drives the Uniformization method
 * for simulating end-conditioned realizations
 * <p/>
 * This work is supported by NSF grant 0856099
 * <p/>
 * Hobolth A and Stone E (2009) Simulation from endpoint-conditioned, continuous-time Markov chains on a finite
 * state space, with applications to molecular evolution. Annals of Applied Statistics, 3, 1204-1231.
 *
 * @author Marc A. Suchard
 */

public class SubordinatedProcess {

    public SubordinatedProcess(double[] Q, int stateCount) {
        this.stateCount = stateCount;
        poissonRate = getMaxRate(Q, stateCount);
        dtmcCache = new ArrayList<double[]>();
        dtmcCache.add(makeIndentityMatrx(stateCount));
        dtmcCache.add(constructDtmcMatrix(Q, stateCount));
        tmp = new double[stateCount];
        this.Q = Q;
    }

    public double getPoissonRate() {
        return poissonRate;
    }

    private double getCachedExp(double x) {
        if (x != cachedXForExp) {
            cachedXForExp = x;
            cachedExpValue = Math.exp(x);
        }
        return cachedExpValue;
    }

    /**
     * Compute the n-step discrete-time transition probabilities
     *
     * @param nSteps which step
     * @return a pointer to the cached matrix
     */

    public double[] getDtmcProbabilities(int nSteps) {
        if (nSteps > dtmcCache.size() - 1) {
            double[] dtmcOneStep = dtmcCache.get(1);
            for (int step = dtmcCache.size() - 1; step <= nSteps; step++) {
                double[] lastDtmcMatrix = dtmcCache.get(step);
                double[] nextDtmcMatrix = new double[stateCount * stateCount];
                MarkovJumpsCore.matrixMultiply(lastDtmcMatrix, dtmcOneStep, stateCount, nextDtmcMatrix);
                dtmcCache.add(nextDtmcMatrix);
            }
        }
        return dtmcCache.get(nSteps);
    }

    /**
     * Find max_i -Q_{ii}
     *
     * @param Q          ctmc rate matrix
     * @param stateCount dim
     * @return max rate
     */

    private double getMaxRate(double[] Q, int stateCount) {
        double max = -Q[0];
        for (int i = 1; i < stateCount; i++) {
            double nextRate = -Q[i * stateCount + i];
            if (nextRate > max) {
                max = nextRate;
            }
        }
        return max;
    }

    /**
     * R =  I + 1/maxRate Q
     *
     * @param lambda     Q
     * @param stateCount dim
     * @return R
     */

    private double[] constructDtmcMatrix(double[] lambda, int stateCount) {
        double[] R = new double[stateCount * stateCount];
        double maxRate = getMaxRate(lambda, stateCount);
        int index = 0;
        for (int i = 0; i < stateCount; i++) {
            for (int j = 0; j < stateCount; j++) {
                R[index] = lambda[index] / maxRate;

                if (i == j) {
                    R[index] += 1;
                }
                index++;
            }
        }
        return R;
    }

    /**
     * Simulate transition times, uniformly distributed before sorting
     *
     * @param timeDuration         T
     * @param totalNumberOfChanges total number of changes
     * @return the transition times of the subordinated process
     */
    public double[] drawTransitionTimes(double timeDuration, int totalNumberOfChanges) {
        double[] times = new double[totalNumberOfChanges];
        for (int i = 0; i < totalNumberOfChanges; i++) {
            times[i] = timeDuration * MathUtils.nextDouble();
        }
        if (times.length > 1) {
            Arrays.sort(times);
        }
        return times;
    }

    /**
     * Simulate the next transition in the subordinated process, equation in remark 7
     *
     * @param currentState         current state of the subordinated process
     * @param endingState          ending state of CTMC
     * @param totalNumberOfChanges number of subordinated changes
     * @param thisChangeNumber     this transition number
     * @return the next state of the subordinated process
     */
    public int drawNextChainState(int currentState, int endingState, int totalNumberOfChanges, int thisChangeNumber) {
        computePdfNextChainState(currentState, endingState, totalNumberOfChanges, thisChangeNumber, tmp);
        return MathUtils.randomChoicePDF(tmp);
    }

    public void computePdfNextChainState(int currentState, int endingState, int totalNumberOfChanges, int thisChangeNumber,
                                         double[] pdf) {
        double[] R = getDtmcProbabilities(1);
        double[] RnMinusI = getDtmcProbabilities(totalNumberOfChanges - thisChangeNumber);

        for (int i = 0; i < stateCount; i++) {
            pdf[i] = R[currentState * stateCount + i] * RnMinusI[i * stateCount + endingState];
//                     / RnMinusIPlus1[currentState * stateCount + endingState] // No need to normalize
        }
    }

    public class Exception extends java.lang.Exception {
        // Nothing special
    }

    /**
     * Simulate the number of transitions in the subordinated process, equation (2.9)
     *
     * @param startingState   starting state of CTMC
     * @param endingState     ending state of CTMC
     * @param time            length of chain
     * @param ctmcProbability the CTMC finite-time transition probability
     * @return the number of transitions in the subordinated process
     * @throws dr.inference.markovjumps.SubordinatedProcess.Exception exception
     */
    public int drawNumberOfChanges(int startingState, int endingState, double time, double ctmcProbability) throws SubordinatedProcess.Exception {
        return drawNumberOfChanges(startingState, endingState, time, ctmcProbability, MathUtils.nextDouble());
    }

    public int drawNumberOfChanges(int startingState, int endingState, double time, double ctmcProbability,
                                   double cutoff) throws SubordinatedProcess.Exception {
        int drawnNumber = -1;
        double cdf = 0;

        double effectiveRate = getPoissonRate() * time;
        double preFactor = getCachedExp(-effectiveRate);
        double scale = 1.0;
        int index = startingState * stateCount + endingState;

        double[] check;
        int maxTries = 1000;
        if (DEBUG) {
            check = new double[maxTries+1];
        }

        while (cutoff >= cdf) {
            drawnNumber++;

            double[] Rn = getDtmcProbabilities(drawnNumber);
            if (drawnNumber > 0) {
                scale *= effectiveRate;
            }
            if (drawnNumber > 1) {
                scale /= (double) drawnNumber;
            }

            cdf += preFactor * scale * Rn[index] / ctmcProbability;

            if (THROW_EXCEPTION) {
                if (drawnNumber == maxTries) {
                    throw new SubordinatedProcess.Exception();
                }
            }

            if (DEBUG) {
                check[drawnNumber] = cdf;
                if (drawnNumber == maxTries) {
                    System.err.println("Start state = " + startingState);
                    System.err.println("End state   = " + endingState);
                    System.err.println("Time        = " + time);
                    System.err.println("CDF         = " + cdf);
                    System.err.println("Cutoff      = " + cutoff);
                    System.err.println("CTMC prob   = " + ctmcProbability);
                    System.err.println("PoissonRate = " + getPoissonRate());

                    double[] distr = computePDFDirectly(startingState, endingState, time, ctmcProbability, drawnNumber);
                    double[] checkCDF = new double[distr.length];
                    double total = 0;
                    for (int i = 0; i < distr.length; i++) {
                        total += distr[i];
                        checkCDF[i] = total;
                    }
                    System.err.println("Direct compute = " + new Vector(distr));
                    System.err.println("Via CDF        = " + new Vector(checkCDF));
                    System.err.println("Check distr    = " + new Vector(check));
                    System.err.println("Q              = " + new Vector(Q));
                    System.err.println("R              = " + new Vector(getDtmcProbabilities(1)));

                    throw new RuntimeException("Likely numerical instability in computing end-conditioned CTMC simulant.");
                }
            }
        }
        return drawnNumber;
    }

    public double[] computePDFDirectly(int startingState, int endingState, double time, double ctmcProbability,
                                       int maxTerm) {
        double[] pdf = new double[maxTerm];

        final double logRateTime = Math.log(getPoissonRate())+ Math.log(time);
        final double logCtmcProbability = Math.log(ctmcProbability);

        for (int n = 0; n < maxTerm; n++) {
            double[] Rn = getDtmcProbabilities(n);
//            pdf[n] = Math.exp(-getPoissonRate() * time) * Math.pow(getPoissonRate() * time, n) /
//                    Math.exp(GammaFunction.lnGamma(n + 1)) * Rn[startingState * stateCount + endingState] /
//                    ctmcProbability;

            pdf[n] = Math.exp(-getPoissonRate() * time + n * logRateTime - GammaFunction.lnGamma(n + 1) +
                    Math.log(Rn[startingState * stateCount + endingState]) - logCtmcProbability);
        }
        return pdf;
    }

    private double[] makeIndentityMatrx(int stateCount) {
        double[] I = new double[stateCount * stateCount];
        for (int i = 0; i < stateCount; i++) {
            I[i * stateCount + i] = 1.0;
        }
        return I;
    }

    private final List<double[]> dtmcCache;
    private final double poissonRate;
    private final int stateCount;
    private final double[] tmp;

    private double cachedXForExp = Double.NaN;
    private double cachedExpValue;

    private static final boolean DEBUG = false;
    private static final boolean THROW_EXCEPTION = true;

    private double[] Q;
}
