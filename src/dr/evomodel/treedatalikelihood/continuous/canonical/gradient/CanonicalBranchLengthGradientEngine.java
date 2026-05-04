package dr.evomodel.treedatalikelihood.continuous.canonical.gradient;

import dr.evomodel.continuous.ou.canonical.CanonicalOUKernel;
import dr.evomodel.continuous.ou.OUProcessModel;

import java.util.Arrays;

public final class CanonicalBranchLengthGradientEngine {

    public void compute(final OUProcessModel processModel,
                 final BranchGradientInputs inputs,
                 final double[] gradT) {
        Arrays.fill(gradT, 0.0);
        final CanonicalOUKernel kernel = processModel.getCanonicalKernel();
        for (int activeIndex = 0; activeIndex < inputs.getActiveBranchCount(); ++activeIndex) {
            gradT[inputs.getActiveChildIndex(activeIndex)] =
                    kernel.contractBranchLengthGradientFlat(
                            inputs.getBranchLength(activeIndex),
                            inputs.getTransitionMatrixAdjointFlat(activeIndex),
                            inputs.getTransitionOffsetAdjoint(activeIndex),
                            inputs.getTransitionCovarianceAdjointFlat(activeIndex),
                            false);
        }
    }
}
