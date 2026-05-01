package dr.evomodel.treedatalikelihood.continuous;

import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalLocalTransitionAdjoints;

final class OUCanonicalBranchGradientScratch {

    final CanonicalLocalTransitionAdjoints localAdjoints;
    final double[] branchMeanGradient;
    final double[] matrixGradient;
    final double[] stationaryMean;
    final double[] compressedNativeGradient;
    final double[][] rotationGradient;
    final double[][] transitionAdjoint;
    final double[][] covarianceAdjoint;
    final double[] parameterRestore;

    OUCanonicalBranchGradientScratch(final int dimension) {
        this.localAdjoints = new CanonicalLocalTransitionAdjoints(dimension);
        this.branchMeanGradient = new double[dimension];
        this.matrixGradient = new double[dimension * dimension];
        this.stationaryMean = new double[dimension];
        this.compressedNativeGradient = new double[dimension + 2 * (dimension / 2)];
        this.rotationGradient = new double[dimension][dimension];
        this.transitionAdjoint = new double[dimension][dimension];
        this.covarianceAdjoint = new double[dimension][dimension];
        this.parameterRestore = new double[Math.max(1, dimension * dimension)];
    }
}
