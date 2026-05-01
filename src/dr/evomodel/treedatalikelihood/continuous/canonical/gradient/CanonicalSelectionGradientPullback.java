package dr.evomodel.treedatalikelihood.continuous.canonical.gradient;

import dr.evomodel.treedatalikelihood.continuous.canonical.*;

import dr.evomodel.treedatalikelihood.continuous.canonical.workspace.*;

import dr.evomodel.continuous.ou.OUProcessModel;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalLocalTransitionAdjoints;

interface CanonicalSelectionGradientPullback {

    void initialize(BranchGradientWorkspace workspace,
                    double[] gradA,
                    double[] gradMu);

    void clearWorkerBuffers(BranchGradientWorkspace workspace,
                            int gradALength,
                            int gradMuLength);

    public void prepareWorkspace(BranchGradientWorkspace workspace);

    void accumulateForBranch(OUProcessModel processModel,
                             BranchGradientInputs inputs,
                             int activeIndex,
                             CanonicalLocalTransitionAdjoints localAdjoints,
                             BranchGradientWorkspace workspace,
                             double[] gradA,
                             double[] gradQ,
                             double[] gradMu);

    void reduceWorker(BranchGradientWorkspace worker,
                      BranchGradientWorkspace reductionWorkspace,
                      double[] gradA);

    void finish(BranchGradientInputs inputs,
                BranchGradientWorkspace workspace,
                double[] gradA);
}
