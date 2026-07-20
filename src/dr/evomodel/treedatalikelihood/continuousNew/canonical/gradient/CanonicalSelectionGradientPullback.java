package dr.evomodel.treedatalikelihood.continuousNew.canonical.gradient;

import dr.evomodel.continuous.ou.OUProcessModel;
import dr.evomodel.treedatalikelihood.continuousNew.canonical.message.CanonicalLocalTransitionAdjoints;
import dr.evomodel.continuous.canonical.workspace.BranchGradientWorkspace;

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
                double[] gradA,
                double[] gradQ,
                double[] gradMu);
}
