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
public class NewTipFullConditionalDistributionDelegate extends
        AbstractFullConditionalDistributionDelegate {

    private static final String NAME_PREFIX = "nFcd";

    public NewTipFullConditionalDistributionDelegate(String name, Tree tree,
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

        TreeTrait<List<NormalSufficientStatistics>> baseTrait =
                new TreeTrait<List<NormalSufficientStatistics>>() {

            public String getTraitName() {
                return getName(name);
            }

            public Intent getIntent() {
                return Intent.NODE;
            }

            public Class getTraitClass() { return List.class; }

            public List<NormalSufficientStatistics> getTrait(Tree t, NodeRef node) {
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

    private static String formatted(List<NormalSufficientStatistics> statistics) {

        StringBuilder sb = new StringBuilder();
        for (NormalSufficientStatistics stat : statistics) {
            sb.append(stat.toVectorizedString()).append(";\t");
        }

        return sb.toString();
    }

    protected List<NormalSufficientStatistics> getTraitForNode(NodeRef node) {

        assert simulationProcess != null;
        assert dimPartial > 0;
        assert dimTrait > 0;

        simulationProcess.cacheSimulatedTraits(node);

        final int numberOfNodes = (node == null) ? tree.getNodeCount() : 1;
        double[] partial = new double[dimPartial * numTraits * numberOfNodes * 4];

        final int index = (node == null) ? -1 : likelihoodDelegate.getActiveNodeIndex(node.getNumber());
        cdi.getPreOrderPartial(index, partial);

        List<NormalSufficientStatistics> statistics = new ArrayList<NormalSufficientStatistics>();

        for (int n = 0; n < numberOfNodes; ++n) {
            for (int i = 0; i < numTraits; ++i) {
                statistics.add(new NormalSufficientStatistics(partial, n * dimPartial, dimTrait,
                        Pd, likelihoodDelegate.getPrecisionType()));
            }
        }

        return statistics;
    }
}
