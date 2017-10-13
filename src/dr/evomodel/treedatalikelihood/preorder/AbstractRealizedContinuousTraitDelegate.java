package dr.evomodel.treedatalikelihood.preorder;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evomodel.continuous.MultivariateDiffusionModel;
import dr.evomodel.treedatalikelihood.ProcessOnTreeDelegate;
import dr.evomodel.treedatalikelihood.continuous.ConjugateRootTraitPrior;
import dr.evomodel.treedatalikelihood.continuous.ContinuousDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.continuous.ContinuousRateTransformation;
import dr.evomodel.treedatalikelihood.continuous.ContinuousTraitPartialsProvider;
import dr.evomodel.treedatalikelihood.continuous.cdi.ContinuousDiffusionIntegrator;

import java.util.Arrays;
import java.util.List;

/**
 * @author Marc A. Suchard
 */
public abstract class AbstractRealizedContinuousTraitDelegate extends ProcessSimulationDelegate.AbstractContinuousTraitDelegate {

    AbstractRealizedContinuousTraitDelegate(String name,
                                            Tree tree,
                                            MultivariateDiffusionModel diffusionModel,
                                            ContinuousTraitPartialsProvider dataModel,
                                            ConjugateRootTraitPrior rootPrior,
                                            ContinuousRateTransformation rateTransformation,
                                            ContinuousDataLikelihoodDelegate likelihoodDelegate) {
        super(name, tree, diffusionModel, dataModel, rootPrior, rateTransformation, likelihoodDelegate);

        sample = new double[dimNode * tree.getNodeCount()];
        tmpEpsilon = new double[dimTrait];
        tmpDrift = new double[dimTrait];
    }

    @Override
    protected void constructTraits(final Helper treeTraitHelper) {

        TreeTrait.DA baseTrait = new TreeTrait.DA() {

            public String getTraitName() {
                return name;
            }

            public Intent getIntent() {
                return Intent.NODE;
            }

            public double[] getTrait(Tree t, NodeRef node) {

                if (t != tree) {  // TODO Write a wrapper class around t if TransformableTree
                    if (t == baseTree) {
                        node = getBaseNode(t, node);
                    } else {
                        throw new RuntimeException("Tree '" + t.getId() + "' and likelihood '" + tree.getId() + "' mismatch");
                    }
                }

                return getTraitForNode(node);
            }
        };

        treeTraitHelper.addTrait(baseTrait);

        TreeTrait.DA tipTrait = new TreeTrait.DA() {

            public String getTraitName() {
                return getTipTraitName(name);
            }

            public Intent getIntent() {
                return Intent.WHOLE_TREE;
            }

            public double[] getTrait(Tree t, NodeRef node) {

                assert t == tree;
                return getTraitForAllTips();
            }
        };

        treeTraitHelper.addTrait(tipTrait);

        TreeTrait.DA tipPrecision = new TreeTrait.DA() {

            public String getTraitName() {
                return getTipPrecisionName(name);
            }

            public Intent getIntent() {
                return Intent.WHOLE_TREE;
            }

            public double[] getTrait(Tree t, NodeRef node) {

                assert t == tree;
                return getPrecisionForAllTips();
            }
        };

        treeTraitHelper.addTrait(tipPrecision);
    }

    public static String getTipTraitName(String name) {
        return "tip." + name;
    }

    private static String getTipPrecisionName(String name) {
        return "precision." + name;
    }

    private double[] getTraitForAllTips() {

        assert simulationProcess != null;

        simulationProcess.cacheSimulatedTraits(null);

        final int length = dimNode * tree.getExternalNodeCount();
        double[] trait = new double[length];
        System.arraycopy(sample, 0, trait, 0, length);

        return trait;
    }

    private double[] getPrecisionForAllTips() {

        assert simulationProcess != null;

        simulationProcess.cacheSimulatedTraits(null);

        final int length = tree.getExternalNodeCount();
        double[] precision = new double[length];

        Arrays.fill(precision, Double.POSITIVE_INFINITY); // TODO

        return precision;
    }

    private double[] getTraitForNode(final NodeRef node) {

        assert simulationProcess != null;

        simulationProcess.cacheSimulatedTraits(null);

        if (node == null) {
            return getTraitForAllTips();
        } else {
            double[] trait = new double[dimNode];
            System.arraycopy(sample, node.getNumber() * dimNode, trait, 0, dimNode);

            return trait;
        }
    }

    public int vectorizeNodeOperations(final List<NodeOperation> nodeOperations, final int[] operations) {

        int k = 0;
        for (ProcessOnTreeDelegate.NodeOperation op : nodeOperations) {

            operations[k    ] = op.getNodeNumber(); // Parent sample
            operations[k + 1] = op.getLeftChild();  // Node sample
            operations[k + 2] = likelihoodDelegate.getActiveNodeIndex(op.getLeftChild());   // Node post-order partial
            operations[k + 3] = likelihoodDelegate.getActiveMatrixIndex(op.getLeftChild()); // Node branch info
            operations[k + 4] = op.getLeftChild() < tree.getExternalNodeCount() ? 1 : 0;    // Is node external?

            k += ContinuousDiffusionIntegrator.OPERATION_TUPLE_SIZE;
        }

        return nodeOperations.size();
    }

    final double[] sample;
    final double[] tmpEpsilon;
    final double[] tmpDrift;
}
