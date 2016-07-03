/*
 * StateHistoryTest.java
 *
 * Copyright (c) 2002-2014 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package test.dr.evomodel.substmodel;

import dr.evomodel.substmodel.nucleotide.HKY;
import dr.evolution.datatype.Nucleotides;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.MarkovJumpsSubstitutionModel;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evomodel.substmodel.UniformizedSubstitutionModel;
import dr.inference.markovjumps.MarkovJumpsCore;
import dr.inference.markovjumps.StateHistory;
import dr.math.LogTricks;
import dr.math.MathUtils;
import dr.math.matrixAlgebra.Vector;
import test.dr.math.MathTestCase;

/**
 * @author Marc A. Suchard
 */
public class StateHistoryTest extends MathTestCase {

    public static final int N = 1000000;

    public void setUp() {

        MathUtils.setSeed(666);

        freqModel = new FrequencyModel(Nucleotides.INSTANCE,
                new double[]{0.45, 0.25, 0.05, 0.25});
        baseModel = new HKY(2.0, freqModel);
        stateCount = baseModel.getDataType().getStateCount();

        lambda = new double[stateCount * stateCount];
        baseModel.getInfinitesimalMatrix(lambda);
        System.out.println("lambda = " + new Vector(lambda));

        markovjumps = new MarkovJumpsSubstitutionModel(baseModel);
    }


    public void testGetLogLikelihood() {
        System.out.println("Start of getLogLikelihood test");

        int startingState = 1;
        int endingState = 1;
        double duration = 0.5;
        int iterations = 5000000;
//        int iterations = 1;

        double[] probs = new double[16];
        baseModel.getTransitionProbabilities(duration, probs);
        double trueProb = probs[startingState * 4 + endingState];
        System.out.println("Tru prob = " + trueProb);

        UniformizedSubstitutionModel uSM = new UniformizedSubstitutionModel(baseModel);
        uSM.setSaveCompleteHistory(true);

        double logProb = Double.NEGATIVE_INFINITY;
        double prob = 0.0;
        double condProb = 0.0;
        for (int i = 0; i < iterations; ++i) {
            uSM.computeCondStatMarkovJumps(startingState, endingState, duration, trueProb);
            StateHistory history = uSM.getStateHistory();
//            System.out.println(history.getEndingTime() - history.getStartingTime());
            assertEquals(history.getEndingTime() - history.getStartingTime(), duration, 10E-3);
            assertEquals(history.getStartingState(), startingState);
            assertEquals(history.getEndingState(), endingState);
            double logLikelihood = history.getLogLikelihood(lambda, 4);
            prob += Math.exp(logLikelihood);
            logProb = LogTricks.logSum(logProb, logLikelihood);
            condProb += Math.exp(-logLikelihood);
//            System.err.println(logLikelihood);
        }
        logProb = Math.exp(logProb - Math.log(iterations));
        prob /= iterations;
        condProb /= iterations;
        System.out.println("Sim prob = " + prob);
        System.out.println("Inv prob = " + (1.0 / condProb));
//        System.out.println("log prob = " + logProb);
        //   System.exit(-1);


        System.out.println();
        System.out.println();

        // Try using unconditioned simulation
        double marginalProb = 0.0;
        double mcProb = 0.0;
        double invMcProb = 0.0;
        int totalTries = 0;
        int i = 0;
        while (i < iterations) {
            startingState = MathUtils.randomChoicePDF(freqModel.getFrequencies());
            StateHistory history = StateHistory.simulateUnconditionalOnEndingState(0, startingState, duration, lambda, 4);
            if (//history.getEndingState() == endingState &&
//                    history.getNumberOfJumps() == randomChoice1
                    true
                    ) {
                marginalProb += 1.0;
                i++;
                double logLike = history.getLogLikelihood(lambda, 4);
                mcProb += Math.exp(logLike);
                invMcProb += Math.exp(-logLike);

//                if (i % 100000 == 0) System.out.println(i);
            }
            totalTries++;
        }
        marginalProb /= totalTries;
        mcProb /= iterations;
        invMcProb /= iterations;
        System.out.println("Sim uncd = " + marginalProb);
        System.out.println("mc  prob = " + mcProb);
        System.out.println("m2  prob = " + (1.0 / invMcProb));

        assertEquals(prob, trueProb);
    }

    public void testFreqDistribution() {

        System.out.println("Start of FreqDistribution test");
        int startingState = 0;
        double duration = 10; // 10 expected substitutions is close to \infty

        double[] freq = new double[stateCount];

        for (int i = 0; i < N; i++) {
            StateHistory simultant = StateHistory.simulateUnconditionalOnEndingState(0.0, startingState, duration,
                    lambda, stateCount);
            freq[simultant.getEndingState()]++;
        }

        for (int i = 0; i < stateCount; i++) {
            freq[i] /= N;
        }

        System.out.println("freq = " + new Vector(freq));
        assertEquals(freq, freqModel.getFrequencies(), 1E-3);
        System.out.println("End of FreqDistribution test\n");
    }

    public void testCounts() {

        System.out.println("State of Counts test");
        int startingState = 2;
        double duration = 0.5;

        int[] counts = new int[stateCount * stateCount];
        double[] expectedCounts = new double[stateCount * stateCount];

        for (int i = 0; i < N; i++) {
            StateHistory simultant = StateHistory.simulateUnconditionalOnEndingState(0.0, startingState, duration,
                    lambda, stateCount);
            simultant.accumulateSufficientStatistics(counts, null);
        }

        for (int i = 0; i < stateCount * stateCount; i++) {
            expectedCounts[i] = (double) counts[i] / (double) N;
        }

        double[] r = new double[stateCount * stateCount];
        double[] joint = new double[stateCount * stateCount];
        double[] analytic = new double[stateCount * stateCount];

        for (int from = 0; from < stateCount; from++) {
            for (int to = 0; to < stateCount; to++) {
                double marginal = 0;
                if (from != to) {
                    MarkovJumpsCore.fillRegistrationMatrix(r, from, to, stateCount);
                    markovjumps.setRegistration(r);
                    markovjumps.computeJointStatMarkovJumps(duration, joint);

                    for (int j = 0; j < stateCount; j++) {
                        marginal += joint[startingState * stateCount + j]; // Marginalize out ending state
                    }
                }
                analytic[from * stateCount + to] = marginal;
            }
        }

        System.out.println("unconditional expected counts = " + new Vector(expectedCounts));
        System.out.println("analytic               counts = " + new Vector(analytic));

        assertEquals(expectedCounts, analytic, 1E-3);
        System.out.println("End of Counts test\n");
    }

    double[] lambda;
    FrequencyModel freqModel;
    SubstitutionModel baseModel;
    MarkovJumpsSubstitutionModel markovjumps;
    int stateCount;


}
