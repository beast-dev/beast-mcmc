package dr.evomodel.treedatalikelihood.preorder;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evomodel.continuous.MultivariateDiffusionModel;
import dr.evomodel.treedatalikelihood.continuous.ConjugateRootTraitPrior;
import dr.evomodel.treedatalikelihood.continuous.ContinuousDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.continuous.ContinuousRateTransformation;
import dr.evomodel.treedatalikelihood.continuous.ContinuousTraitPartialsProvider;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Marc A. Suchard
 */
public class BranchConditionalDistributionDelegate extends
        AbstractFullConditionalDistributionDelegate {

    private static final String NAME_PREFIX = "bcd";

    public BranchConditionalDistributionDelegate(String name, Tree tree,
                                                 MultivariateDiffusionModel diffusionModel,
                                                 ContinuousTraitPartialsProvider dataModel,
                                                 ConjugateRootTraitPrior rootPrior,
                                                 ContinuousRateTransformation rateTransformation,
                                                 ContinuousDataLikelihoodDelegate likelihoodDelegate) {
        super(name, tree, diffusionModel, dataModel, rootPrior, rateTransformation, likelihoodDelegate);
    }

    public static String getName(String name) {
        return NAME_PREFIX + "." + name;
    }

    @Override
    protected void constructTraits(Helper treeTraitHelper) {

        TreeTrait<List<BranchSufficientStatistics>> baseTrait =
                new TreeTrait<List<BranchSufficientStatistics>>() {

            public String getTraitName() {
                return getName(name);
            }

            public Intent getIntent() {
                return Intent.BRANCH;
            }

            public Class getTraitClass() { return List.class; }

            public List<BranchSufficientStatistics> getTrait(Tree t, NodeRef node) {
                assert (tree == t);

                return getTraitForNode(node);
            }

            public String getTraitString(Tree tree, NodeRef node) {
                return formatted(getTrait(tree, node));
            }

            public boolean getLoggable() {
                return isLoggable();
            }
        };

        treeTraitHelper.addTrait(baseTrait);
    }

    private static String formatted(List<BranchSufficientStatistics> statistics) {

        StringBuilder sb = new StringBuilder();
        for (BranchSufficientStatistics stat : statistics) {
            sb.append(stat.toVectorizedString()).append(";\t");
        }

        return sb.toString();
    }

    protected List<BranchSufficientStatistics> getTraitForNode(NodeRef node) {

        assert simulationProcess != null;
        assert dimPartial > 0;
        assert dimTrait > 0;

        simulationProcess.cacheSimulatedTraits(node);

        int numberOfNodes = (node == null) ? likelihoodDelegate.getPartialBufferCount() : 1;
        double[] childPartial = new double[dimPartial * numTraits * numberOfNodes];
        double[] parentPartial = new double[dimPartial * numTraits * numberOfNodes];
        double[] branchPrecision = new double[dimTrait * dimTrait * numTraits * numberOfNodes];
        double[] branchDisplacement = new double[dimTrait * numberOfNodes];

        int nodeNumber = (node == null) ? -1 : likelihoodDelegate.getActiveNodeIndex(node.getNumber());
        int branchNumber = (node == null) ? -1 : likelihoodDelegate.getActiveMatrixIndex(node.getNumber());

        List<BranchSufficientStatistics> statistics = new ArrayList<BranchSufficientStatistics>();

        cdi.getPostOrderPartial(nodeNumber, childPartial);
        cdi.getBranchMatrices(branchNumber, branchPrecision, branchDisplacement);
        cdi.getPreOrderPartial(nodeNumber, parentPartial);


        
        if (node == null) {
            for (int n = 0; n < tree.getNodeCount(); ++n) {

                NodeRef tNode = tree.getNode(n);
                int nodeIndex = likelihoodDelegate.getActiveNodeIndex(tNode.getNumber());
                int branchIndex = likelihoodDelegate.getActiveMatrixIndex(tNode.getNumber());
                addOneNode(statistics, childPartial, parentPartial, branchPrecision, branchDisplacement,
                        nodeIndex, branchIndex);
            }
        } else {
            addOneNode(statistics, childPartial, parentPartial, branchPrecision, branchDisplacement, 0, 0);
        }

        return statistics;
    }

    private void addOneNode(List<BranchSufficientStatistics> statistics,
                            double[] childPartial, double[] parentPartial,
                            double[] branchPrecision, double[] branchDisplacement,
                            int nodeIndex, int branchIndex) {
        if (numTraits > 1) throw new IllegalArgumentException("Not yet implemented");

//        for (int i = 0; i < numTraits; ++i) {  // TODO Enable for > 1 numTraits
            statistics.add(
                    new BranchSufficientStatistics(
                            new NormalSufficientStatistics(childPartial, nodeIndex, dimTrait,
                                    Pd, likelihoodDelegate.getPrecisionType()),
                            new NormalSufficientStatistics(branchDisplacement, branchPrecision, branchIndex, dimTrait,
                                    Pd, likelihoodDelegate.getPrecisionType()),
                            new NormalSufficientStatistics(parentPartial, nodeIndex, dimTrait,
                                    Pd, likelihoodDelegate.getPrecisionType())
                    ));
//            offset +=  dimPartial;
//        }
    }
}
