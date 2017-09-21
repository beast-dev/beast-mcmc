package dr.evomodel.treedatalikelihood.preorder;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.continuous.MultivariateDiffusionModel;
import dr.evomodel.treedatalikelihood.continuous.ConjugateRootTraitPrior;
import dr.evomodel.treedatalikelihood.continuous.ContinuousDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.continuous.ContinuousRateTransformation;
import dr.evomodel.treedatalikelihood.continuous.ContinuousTraitPartialsProvider;

import java.util.Arrays;

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
                                            BranchRateModel rateModel,
                                            ContinuousDataLikelihoodDelegate likelihoodDelegate) {
        super(name, tree, diffusionModel, dataModel, rootPrior, rateTransformation, rateModel, likelihoodDelegate);

        sample = new double[dimNode * tree.getNodeCount()];
        tmpEpsilon = new double[dimTrait];
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

    public static String getTipPrecisionName(String name) {
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

    protected final double[] sample;
    protected final double[] tmpEpsilon;
}
