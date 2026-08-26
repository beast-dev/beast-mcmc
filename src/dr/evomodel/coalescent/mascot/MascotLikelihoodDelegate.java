/*
 * MascotLikelihoodDelegate.java
 *
 * Copyright (c) 2002-2026 Alexei Drummond, Andrew Rambaut and Marc Suchard
 */

package dr.evomodel.coalescent.mascot;

/**
 * MASCOT-native likelihood delegate boundary.
 */
public interface MascotLikelihoodDelegate {

    double calculateLikelihood(MascotPreparedInput prepared, double[] theta, double[] branchRates,
                               boolean checkProbabilities);

    double calculateLikelihoodAndDerivatives(MascotPreparedInput prepared, double[] theta, double[] branchRates,
                                             double[] gradientOut, double[] clockGradientOut,
                                             double[] ancestralStateScoresOut, boolean checkProbabilities);

    final class Result {
        public final double logLikelihood;
        public final double[] gradient;
        public final double[] clockGradient;
        public final double[] rootProbabilities;
        public final int[] activeLineages;
        public final double[] ancestralStateScores;

        Result(double logLikelihood, double[] gradient, double[] clockGradient,
               double[] rootProbabilities, int[] activeLineages, double[] ancestralStateScores) {
            this.logLikelihood = logLikelihood;
            this.gradient = gradient;
            this.clockGradient = clockGradient;
            this.rootProbabilities = rootProbabilities;
            this.activeLineages = activeLineages;
            this.ancestralStateScores = ancestralStateScores;
        }
    }

    final class NumericalException extends RuntimeException {
        public NumericalException(String message) {
            super(message);
        }
    }
}
