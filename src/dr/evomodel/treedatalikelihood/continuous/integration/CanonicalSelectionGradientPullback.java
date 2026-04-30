package dr.evomodel.treedatalikelihood.continuous.integration;

import dr.evomodel.continuous.ou.OUProcessModel;

interface CanonicalSelectionGradientPullback {

    void initialize(BranchGradientWorkspace workspace,
                    double[] gradA,
                    double[] gradMu);

    void clearWorkerBuffers(BranchGradientWorkspace workspace,
                            int gradALength,
                            int gradMuLength);

    void prepareWorkspace(BranchGradientWorkspace workspace);

    void accumulateForBranch(OUProcessModel processModel,
                             BranchGradientInputs inputs,
                             int activeIndex,
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
