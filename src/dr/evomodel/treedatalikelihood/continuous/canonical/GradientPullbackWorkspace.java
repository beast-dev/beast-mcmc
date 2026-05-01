package dr.evomodel.treedatalikelihood.continuous.canonical;

import dr.evomodel.continuous.ou.canonical.CanonicalBranchWorkspace;
import dr.evomodel.continuous.ou.canonical.CanonicalPreparedTransitionCapability;
import dr.evomodel.continuous.ou.orthogonalblockdiagonal.OrthogonalBlockBranchGradientWorkspace;
import dr.evomodel.continuous.ou.orthogonalblockdiagonal.OrthogonalBlockCanonicalParameterization;

import java.util.Arrays;

final class GradientPullbackWorkspace {
    final double[] orthogonalStationaryMeanScratch;
    final double[] orthogonalCompressedGradientScratch;
    final double[] orthogonalNativeGradientScratch;
    final double[] orthogonalRotationGradientFlatScratch;
    final double[][] orthogonalRotationGradientScratch;
    final double[][] covariance2;
    final double[] transitionMatrixFlat;
    final double[] covarianceAdjointFlat;
    final double[] matrixProductFlat;
    final double[] localGradientA;
    final double[] localGradientQ;
    final double[] localGradientMuVector;
    final double[] localGradientMuScalar;
    CanonicalBranchWorkspace specializedBranchWorkspace;

    GradientPullbackWorkspace(final int dim) {
        this.orthogonalStationaryMeanScratch = new double[dim];
        this.orthogonalCompressedGradientScratch = new double[dim + 2 * (dim / 2)];
        this.orthogonalNativeGradientScratch = new double[((dim & 1) == 1 ? 1 : 0) + 3 * (dim / 2)];
        this.orthogonalRotationGradientFlatScratch = new double[dim * dim];
        this.orthogonalRotationGradientScratch = new double[dim][dim];
        this.covariance2 = new double[dim][dim];
        this.transitionMatrixFlat = new double[dim * dim];
        this.covarianceAdjointFlat = new double[dim * dim];
        this.matrixProductFlat = new double[dim * dim];
        this.localGradientA = new double[dim * dim];
        this.localGradientQ = new double[dim * dim];
        this.localGradientMuVector = new double[dim];
        this.localGradientMuScalar = new double[1];
    }

    double[] localGradientMu(final int gradientLength, final int dim) {
        if (gradientLength == 1) {
            return localGradientMuScalar;
        }
        if (gradientLength == dim) {
            return localGradientMuVector;
        }
        throw new IllegalArgumentException(
                "Stationary-mean gradient length must be 1 or " + dim + ", found " + gradientLength);
    }

    void clearLocalGradientBuffers(final int gradALength,
                                   final int gradMuLength,
                                   final int dim,
                                   final boolean orthogonalSelection,
                                   final int compressedGradientLength) {
        Arrays.fill(localGradientA, 0, gradALength, 0.0);
        Arrays.fill(localGradientQ, 0.0);
        Arrays.fill(localGradientMu(gradMuLength, dim), 0.0);
        if (orthogonalSelection) {
            Arrays.fill(orthogonalCompressedGradientScratch, 0, compressedGradientLength, 0.0);
            Arrays.fill(orthogonalRotationGradientFlatScratch, 0.0);
            for (double[] row : orthogonalRotationGradientScratch) {
                Arrays.fill(row, 0.0);
            }
        }
    }

    OrthogonalBlockBranchGradientWorkspace
    ensureOrthogonalBranchWorkspace(final OrthogonalBlockCanonicalParameterization orthogonalSelection) {
        if (specializedBranchWorkspace == null) {
            specializedBranchWorkspace = orthogonalSelection.createBranchGradientWorkspace();
        }
        return (OrthogonalBlockBranchGradientWorkspace) specializedBranchWorkspace;
    }

    CanonicalBranchWorkspace
    ensureSpecializedBranchWorkspace(final CanonicalPreparedTransitionCapability selection) {
        if (specializedBranchWorkspace == null) {
            specializedBranchWorkspace = selection.createBranchWorkspace();
        }
        return specializedBranchWorkspace;
    }
}
