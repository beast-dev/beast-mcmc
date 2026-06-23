package dr.evomodel.treedatalikelihood.continuous.canonical.gradient;

import dr.evomodel.continuous.ou.canonical.CanonicalOUKernel;
import dr.evomodel.continuous.ou.OUProcessModel;
import dr.evomodel.treedatalikelihood.continuous.canonical.CanonicalOUTransitionProvider;

import java.util.Arrays;

public final class CanonicalBranchLengthGradientEngine {

    public void compute(final CanonicalOUTransitionProvider provider,
                        final BranchGradientInputs inputs,
                        final double[] gradT) {
        Arrays.fill(gradT, 0.0);
        final OUProcessModel homogeneousProcessModel =
                provider.getProcessModelCount() == 1 ? provider.getProcessModel(0) : null;
        for (int activeIndex = 0; activeIndex < inputs.getActiveBranchCount(); ++activeIndex) {
            final OUProcessModel processModel = homogeneousProcessModel == null
                    ? provider.getProcessModel(inputs.getProcessIndex(activeIndex))
                    : homogeneousProcessModel;
            final CanonicalOUKernel kernel = processModel.getCanonicalKernel();
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
