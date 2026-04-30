package dr.evomodel.treedatalikelihood.continuous.integration;

import dr.evomodel.continuous.ou.OUProcessModel;

import java.util.Arrays;

final class CanonicalBranchLengthGradientEngine {

    void compute(final OUProcessModel processModel,
                 final BranchGradientInputs inputs,
                 final double[] gradT) {
        Arrays.fill(gradT, 0.0);
        for (int activeIndex = 0; activeIndex < inputs.getActiveBranchCount(); ++activeIndex) {
            gradT[inputs.getActiveChildIndex(activeIndex)] =
                    processModel.contractBranchLengthGradientFlat(
                            inputs.getBranchLength(activeIndex),
                            inputs.getLocalAdjoints(activeIndex).dLogL_dF,
                            inputs.getLocalAdjoints(activeIndex).dLogL_df,
                            inputs.getLocalAdjoints(activeIndex).dLogL_dOmega,
                            false);
        }
    }
}
