package dr.evomodel.treedatalikelihood.markovjumps;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evolution.tree.TreeTraitProvider;
import dr.evomodel.treedatalikelihood.preorder.AbstractRealizedDiscreteTraitDelegate;
import dr.evomodel.treelikelihood.MarkovJumpsTraitProvider;
import dr.inference.markovjumps.StateHistory;
import dr.inference.model.Parameter;

public class MarkovJumpRewardAddOn implements AbstractRealizedDiscreteTraitDelegate.RealizedDiscreteAddOn {

    @SuppressWarnings("unused")
    public MarkovJumpRewardAddOn(String name,
                                 Tree tree,
                                 AbstractRealizedDiscreteTraitDelegate ancestralTraitDelegate,
                                 CompleteHistoryAddOn completeHistoryAddOn,
                                 Parameter register,
                                 MarkovJumpsTraitProvider.ValueScaling valueScaling) {

        this.name = name + "." + register.getParameterName();
        this.ancestralTraitDelegate = ancestralTraitDelegate;
        this.stateCount = ancestralTraitDelegate.getStateCount();
        this.completeHistoryAddOn = completeHistoryAddOn;
        this.register = register;

        if (!(register.getDimension() == stateCount * stateCount ||
                register.getDimension() == stateCount)) {
            throw new IllegalArgumentException("Invalid register dimension");
        }

//        this.valueScaling = valueScaling;
        this.expectedJumpReward = new double[2 * tree.getNodeCount() - 2][];
    }

    private final AbstractRealizedDiscreteTraitDelegate ancestralTraitDelegate;
    private final CompleteHistoryAddOn completeHistoryAddOn;
    private final Parameter register;

    private final double[][] expectedJumpReward;

    private final String name;
    private final int stateCount;

//    private final MarkovJumpsTraitProvider.ValueScaling valueScaling;

    @Override
    public void constructTraits(TreeTraitProvider.Helper treeTraitHelper) {

        TreeTrait<double[]> jumpRewardTrait = new TreeTrait.DA() {

            @Override
            public String getTraitName() {
                return name;
            }

            @Override
            public Intent getIntent() {
                return Intent.BRANCH;
            }

            @Override
            public double[] getTrait(Tree tree, NodeRef node) {
                ancestralTraitDelegate.getStatesForNode(tree, node);
                return expectedJumpReward[node.getNumber()];
            }
        };

        treeTraitHelper.addTrait(jumpRewardTrait);
    }

    @Override
    public void hookCalculation(int node, int parent, double[] probabilities) {
        StateHistory[] history = completeHistoryAddOn.getRawHistory(node);

        if (expectedJumpReward[node] == null) {
            expectedJumpReward[node] = new double[history.length];
        }

        double[] r = register.getParameterValues();
        for (int i = 0; i < history.length; ++i) {
            final double value;
            if (r.length == stateCount * stateCount) {
                value = history[i].getTotalRegisteredCounts(r);
            } else {
                value = history[i].getTotalReward(r);
            }

            // TODO Possible valueScaling here (depends on categoryRate and branchRate

            expectedJumpReward[node][i] = value;
        }
    }
}
