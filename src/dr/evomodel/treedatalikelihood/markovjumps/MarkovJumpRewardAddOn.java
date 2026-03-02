package dr.evomodel.treedatalikelihood.markovjumps;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evolution.tree.TreeTraitProvider;
import dr.evomodel.substmodel.MarkovJumpsSubstitutionModel;
import dr.evomodel.treedatalikelihood.preorder.AbstractRealizedDiscreteTraitDelegate;
import dr.evomodel.treelikelihood.MarkovJumpsTraitProvider;
import dr.inference.markovjumps.StateHistory;
import dr.inference.model.Parameter;

public class MarkovJumpRewardAddOn implements AbstractRealizedDiscreteTraitDelegate.RealizedDiscreteAddOn {

    public static final String MARGINAL_RATE_SUFFIX = ".marginal.rate";

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
        this.register = register.getParameterValues();

        if (!(register.getDimension() == stateCount * stateCount ||
                register.getDimension() == stateCount)) {
            throw new IllegalArgumentException("Invalid register dimension");
        }

//        this.valueScaling = valueScaling;
        this.expectedJumpReward = new double[tree.getNodeCount()][];

        this.marginalRate = new double[completeHistoryAddOn.getMarkovJumpsSubstitutionModelCount()];
        this.marginalRateKnown = false;
    }

    private final AbstractRealizedDiscreteTraitDelegate ancestralTraitDelegate;
    private final CompleteHistoryAddOn completeHistoryAddOn;
    private final double[] register;

    private final double[][] expectedJumpReward;

    private final String name;
    private final int stateCount;

    private final double[] marginalRate;
    private boolean marginalRateKnown;
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

        TreeTrait<double[]> unconditionalRateTrait = new TreeTrait.DA() {

            @Override
            public String getTraitName() {
                return name + MARGINAL_RATE_SUFFIX;
            }

            @Override
            public Intent getIntent() {
                return Intent.WHOLE_TREE;
            }

            @Override
            public double[] getTrait(Tree tree, NodeRef node) {
                updateMarginalRate();
                return marginalRate;
            }
        };

        treeTraitHelper.addTrait(unconditionalRateTrait);
    }

    @Override
    public void hookCalculation(int node, int parent, double[] probabilities) {
        StateHistory[] history = completeHistoryAddOn.getRawHistory(node);

        if (expectedJumpReward[node] == null) {
            expectedJumpReward[node] = new double[history.length];
        }

        for (int i = 0; i < history.length; ++i) {
            final double value;
            if (register.length == stateCount * stateCount) {
                value = history[i].getTotalRegisteredCounts(register);
            } else {
                value = history[i].getTotalReward(register);
            }

            // TODO Possible valueScaling here (depends on categoryRate and branchRate

            expectedJumpReward[node][i] = value;
        }
    }

    public void makeDirty() {
        marginalRateKnown = false;
    }

//    private int getMarginalRate(int index) {
    private void updateMarginalRate() {
        if (!marginalRateKnown) {
            for (int i = 0; i < completeHistoryAddOn.getMarkovJumpsSubstitutionModelCount(); ++i) {
                MarkovJumpsSubstitutionModel mj = completeHistoryAddOn.getMarkovJumpsSubstitutionModel(i);
                mj.setRegistration(register);
                marginalRate[i] = mj.getMarginalRate();
                mj.setRegistration(completeHistoryAddOn.getFullRegistration());
            }
            marginalRateKnown = true;
        }
//        return marginalRate[index];
    }
}
