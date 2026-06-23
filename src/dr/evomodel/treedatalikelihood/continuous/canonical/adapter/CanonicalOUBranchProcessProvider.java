package dr.evomodel.treedatalikelihood.continuous.canonical.adapter;

import dr.evomodel.continuous.ou.OUProcessModel;
import dr.evomodel.treedatalikelihood.continuous.ContinuousBranchTimeModel;
import dr.evomodel.treedatalikelihood.continuous.canonical.CanonicalNativeBranchMessageProvider;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianTransition;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;

/**
 * Branch-process boundary for canonical OU tree likelihoods.
 *
 * <p>The current implementation is homogeneous, but this interface keeps tree traversal
 * and transition caching pointed at a branch-process service rather than directly at
 * one hard-coded OU process model. Heterogeneous providers can later map branches or
 * branch segments to different process models while preserving the same tree-side API.
 */
interface CanonicalOUBranchProcessProvider {

    int getDimension();

    int getNodeCount();

    int getProcessModelCount();

    OUProcessModel getProcessModel(int processIndex);

    int getSegmentCount(int childNodeIndex);

    int getSegmentProcessIndex(int childNodeIndex, int segmentIndex);

    double getSegmentEffectiveTime(int childNodeIndex, int segmentIndex);

    OUProcessModel getProcessModel();

    ContinuousBranchTimeModel getBranchTimeModel();

    double getEffectiveBranchLength(int childNodeIndex);

    void fillTraitCovarianceFlat(double[] out);

    void fillCanonicalTransitionForLength(double branchLength, CanonicalGaussianTransition out);

    void fillTransitionMatrixFlat(double branchLength, double[] out);

    void fillTransitionOffset(double branchLength, double[] out);

    void fillTransitionCovarianceFlat(double branchLength, double[] out);

    void fillTransitionMomentsFlat(double branchLength,
                                   double[] transitionMatrixOut,
                                   double[] transitionOffsetOut,
                                   double[] transitionCovarianceOut);

    CanonicalNativeBranchMessageProvider createNativeBranchMessageProvider(int nodeCount);

    void registerWith(AbstractModel owner);

    void refreshSnapshot();

    void beginPhase();

    void endPhase();

    CanonicalTransitionCacheInvalidationReason classifyModelChange(Model model);

    CanonicalTransitionCacheInvalidationReason classifyVariableChange(Variable variable,
                                                                      Parameter.ChangeType type);
}
