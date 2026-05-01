package dr.evomodel.treedatalikelihood.continuous.canonical;

import dr.evomodel.continuous.ou.canonical.CanonicalBranchWorkspace;
import dr.evomodel.continuous.ou.canonical.CanonicalPreparedTransitionCapability;
import dr.evomodel.continuous.ou.orthogonalblockdiagonal.OrthogonalBlockBranchGradientWorkspace;
import dr.evomodel.continuous.ou.orthogonalblockdiagonal.OrthogonalBlockCanonicalParameterization;

import java.util.Arrays;

final class OrthogonalBlockGradientWorkspace {

    final double[] stationaryMeanScratch;
    final double[] compressedGradientScratch;
    final double[] nativeGradientScratch;
    final double[] rotationGradientFlatScratch;
    final double[][] rotationGradientScratch;
    private CanonicalBranchWorkspace specializedBranchWorkspace;

    OrthogonalBlockGradientWorkspace(final int dim) {
        if (dim < 1) throw new IllegalArgumentException("dim must be >= 1");
        this.stationaryMeanScratch = new double[dim];
        this.compressedGradientScratch = new double[dim + 2 * (dim / 2)];
        this.nativeGradientScratch = new double[((dim & 1) == 1 ? 1 : 0) + 3 * (dim / 2)];
        this.rotationGradientFlatScratch = new double[dim * dim];
        this.rotationGradientScratch = new double[dim][dim];
    }

    void clearSpecializedBuffers(final int compressedGradientLength) {
        Arrays.fill(compressedGradientScratch, 0, compressedGradientLength, 0.0);
        for (double[] row : rotationGradientScratch) {
            Arrays.fill(row, 0.0);
        }
        Arrays.fill(rotationGradientFlatScratch, 0.0);
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
