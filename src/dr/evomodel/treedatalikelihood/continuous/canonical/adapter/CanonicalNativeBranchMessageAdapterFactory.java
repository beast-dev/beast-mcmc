package dr.evomodel.treedatalikelihood.continuous.canonical.adapter;

import dr.evomodel.continuous.ou.SelectionMatrixParameterization;
import dr.evomodel.continuous.ou.kernel.DiagonalOUBranchTransformKernel;
import dr.evomodel.continuous.ou.kernel.OUBranchKernel;
import dr.evomodel.continuous.ou.kernel.OUBranchKernelCapability;
import dr.evomodel.treedatalikelihood.continuous.ContinuousBranchTimeModel;
import dr.evomodel.treedatalikelihood.continuous.canonical.CanonicalNativeBranchMessageProvider;
import dr.inference.model.MatrixParameterInterface;
import dr.inference.model.Parameter;

final class CanonicalNativeBranchMessageAdapterFactory {

    private CanonicalNativeBranchMessageAdapterFactory() {
    }

    static CanonicalNativeBranchMessageProvider create(
            final SelectionMatrixParameterization selection,
            final int nodeCount,
            final MatrixParameterInterface diffusionMatrix,
            final Parameter stationaryMean,
            final ContinuousBranchTimeModel branchTimeModel) {
        if (!(selection instanceof OUBranchKernelCapability)) {
            return null;
        }
        final OUBranchKernel kernel =
                ((OUBranchKernelCapability) selection).createBranchKernel();
        final DiagonalOUBranchTransformKernel diagonalKernel =
                kernel.getCapability(DiagonalOUBranchTransformKernel.class);
        if (diagonalKernel != null) {
            return new DiagonalOUCanonicalBranchMessageAdapter(
                    diagonalKernel,
                    nodeCount,
                    diffusionMatrix,
                    stationaryMean,
                    branchTimeModel);
        }
        return null;
    }
}
