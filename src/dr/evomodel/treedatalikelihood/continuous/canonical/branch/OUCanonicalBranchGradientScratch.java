package dr.evomodel.treedatalikelihood.continuous.canonical.branch;

import dr.evomodel.treedatalikelihood.continuous.CanonicalDebugOptions;
import dr.evomodel.treedatalikelihood.continuous.CanonicalGradientFallbackPolicy;
import dr.evomodel.treedatalikelihood.continuous.OUGaussianBranchTransitionProvider;

import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalLocalTransitionAdjoints;

final class OUCanonicalBranchGradientScratch {

    final CanonicalLocalTransitionAdjoints localAdjoints;
    final double[] branchMeanGradient;
    final double[] matrixGradient;
    final double[] transitionSelectionGradient;
    final double[] covarianceSelectionGradient;
    final double[] totalSelectionGradient;
    final double[] selectionComponentGradient;
    final double[] stationaryMean;
    final double[] compressedNativeGradient;
    final double[] rotationGradient;
    final double[] parameterRestore;

    OUCanonicalBranchGradientScratch(final int dimension) {
        this.localAdjoints = new CanonicalLocalTransitionAdjoints(dimension);
        this.branchMeanGradient = new double[dimension];
        this.matrixGradient = new double[dimension * dimension];
        this.transitionSelectionGradient = new double[dimension * dimension];
        this.covarianceSelectionGradient = new double[dimension * dimension];
        this.totalSelectionGradient = new double[dimension * dimension];
        this.selectionComponentGradient = new double[3 * dimension * dimension];
        this.stationaryMean = new double[dimension];
        this.compressedNativeGradient = new double[dimension + 2 * (dimension / 2)];
        this.rotationGradient = new double[dimension * dimension];
        this.parameterRestore = new double[Math.max(1, dimension * dimension)];
    }
}
