package dr.evomodel.treedatalikelihood.continuous.canonical.workspace;

import dr.evomodel.continuous.ou.canonical.CanonicalBranchWorkspace;
import dr.evomodel.continuous.ou.canonical.CanonicalPreparedTransitionCapability;

import java.util.Arrays;

public final class SpecializedGradientWorkspace {

    public final double[] stationaryMeanScratch;
    public final double[] compressedGradientScratch;
    public final double[] nativeGradientScratch;
    public final double[] rotationGradientFlatScratch;
    public final double[] diffusionGradientDBasisScratch;
    public final double[] meanGradientDBasisScratch;
    private CanonicalBranchWorkspace specializedBranchWorkspace;

    public SpecializedGradientWorkspace(final int dim) {
        if (dim < 1) throw new IllegalArgumentException("dim must be >= 1");
        this.stationaryMeanScratch = new double[dim];
        this.compressedGradientScratch = new double[dim + 2 * (dim / 2)];
        this.nativeGradientScratch = new double[((dim & 1) == 1 ? 1 : 0) + 3 * (dim / 2)];
        this.rotationGradientFlatScratch = new double[dim * dim];
        this.diffusionGradientDBasisScratch = new double[dim * dim];
        this.meanGradientDBasisScratch = new double[dim];
    }

    public void clearSpecializedBuffers(final int compressedGradientLength) {
        Arrays.fill(compressedGradientScratch, 0, compressedGradientLength, 0.0);
        Arrays.fill(rotationGradientFlatScratch, 0.0);
        Arrays.fill(diffusionGradientDBasisScratch, 0.0);
        Arrays.fill(meanGradientDBasisScratch, 0.0);
    }

    public CanonicalBranchWorkspace
    ensureSpecializedBranchWorkspace(final CanonicalPreparedTransitionCapability selection) {
        if (specializedBranchWorkspace == null) {
            specializedBranchWorkspace = selection.createBranchWorkspace();
        }
        return specializedBranchWorkspace;
    }
}
