package dr.evomodel.treedatalikelihood.continuous.canonical.gradient;

/**
 * Supplies branch-local pre/post messages for canonical gradient preparation.
 */
public interface BranchGradientMessageProvider {

    void fillBranchMessages(int childNodeIndex, BranchGradientMessageBundle out);

    double getRootDiffusionScale();

    void fillRootMessages(BranchGradientInputs out);
}
