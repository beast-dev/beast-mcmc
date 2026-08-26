/*
 * MascotLikelihoodBackend.java
 *
 * Copyright (c) 2002-2026 Alexei Drummond, Andrew Rambaut and Marc Suchard
 */

package dr.evomodel.coalescent.mascot;

/**
 * Backend boundary for MASCOT likelihood evaluation.
 */
public interface MascotLikelihoodBackend {

    double evaluateLikelihood(MascotPreparedInput prepared, double[] theta, double[] branchRates,
                              boolean checkProbabilities);

    double evaluateInto(MascotPreparedInput prepared, double[] theta, double[] branchRates,
                        double[] gradientOut, double[] clockGradientOut, double[] ancestralStateScoresOut,
                        boolean checkProbabilities);
}
