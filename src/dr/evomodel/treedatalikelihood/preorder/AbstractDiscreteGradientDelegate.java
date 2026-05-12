package dr.evomodel.treedatalikelihood.preorder;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evomodel.siteratemodel.SiteRateModel;
import dr.evomodel.treedatalikelihood.DiscreteDataLikelihoodDelegate;
import dr.inference.model.Model;

import java.util.List;

/**
 * BEAGLE-free analogue of AbstractBeagleGradientDelegate.
 *
 * This delegate does not maintain its own pre-order buffers.
 * Instead, it asks DiscreteDataLikelihoodDelegate to ensure that
 * post-order and pre-order caches are current, and then computes
 * the requested whole-tree gradient/Hessian-like statistics from them.
 */
/*
* @author Filippo Monti
*/
public abstract class AbstractDiscreteGradientDelegate extends ProcessSimulationDelegate.AbstractDelegate {

    private static final String GRADIENT_TRAIT_NAME = "Gradient";
    private static final String HESSIAN_TRAIT_NAME = "Hessian";

    protected final Tree tree;
    protected final DiscreteDataLikelihoodDelegate likelihoodDelegate;
    protected final SiteRateModel siteRateModel;

    protected final int patternCount;
    protected final int stateCount;
    protected final int categoryCount;

    protected double[] gradient;

    protected boolean substitutionProcessKnown;

    private static final boolean COUNT_TOTAL_OPERATIONS = true;
    private long simulateCount = 0;
    private long getTraitCount = 0;

    protected AbstractDiscreteGradientDelegate(String name,
                                               Tree tree,
                                               DiscreteDataLikelihoodDelegate likelihoodDelegate) {
        super(name, tree);

        this.tree = tree;
        this.likelihoodDelegate = likelihoodDelegate;
        this.siteRateModel = likelihoodDelegate.getSiteRateModel();

        this.patternCount = likelihoodDelegate.getPatternCount();
        this.stateCount = likelihoodDelegate.getStateCount();
        this.categoryCount = likelihoodDelegate.getCategoryCount();

        likelihoodDelegate.addModelListener(this);
        likelihoodDelegate.addModelRestoreListener(this);

        this.substitutionProcessKnown = false;
    }

    protected abstract int getGradientLength();

    protected abstract void getNodeDerivatives(Tree tree, double[] first, double[] second);

    protected String getGradientTraitName() {
        return GRADIENT_TRAIT_NAME;
    }

    protected String getHessianTraitName() {
        return HESSIAN_TRAIT_NAME;
    }

    protected double[] getGradient(NodeRef node) {
        if (COUNT_TOTAL_OPERATIONS) {
            ++getTraitCount;
        }

        simulationProcess.cacheSimulatedTraits(node);
        return gradient.clone();
    }

    protected double[] getHessian(Tree tree, NodeRef node) {
        simulationProcess.cacheSimulatedTraits(node);

        double[] second = new double[getGradientLength()];
        getNodeDerivatives(tree, null, second);

        return second;
    }

    @Override
    protected void setupStatistics() {
        // Make sure the discrete delegate has current post-order / pre-order caches.
        likelihoodDelegate.ensurePreOrderComputed();

        if (gradient == null) {
            gradient = new double[getGradientLength()];
        }

        getNodeDerivatives(tree, gradient, null);
    }

    @Override
    protected void simulateRoot(int rootNumber) {
        // Nothing to do.
        // All required pre/post-order state lives inside DiscreteDataLikelihoodDelegate.
    }

    @Override
    protected void simulateNode(int v0, int v1, int v2, int v3, int v4) {
        // Nothing to do.
        // We do not simulate node-by-node here; setupStatistics() already triggered
        // the required cache preparation and gradient computation.
    }

    @Override
    public int vectorizeNodeOperations(List<NodeOperation> nodeOperations, int rootNodeNumber, int[] operations) {
        // No explicit traversal operations are needed here because the discrete
        // delegate computes and caches pre/post-order internally.
        return 0;
    }

    @Override
    public int getSingleOperationSize() {
        // Unused because vectorizeNodeOperations(...) returns 0.
        return 1;
    }

    @Override
    public void modelChangedEvent(Model model, Object object, int index) {
        substitutionProcessKnown = false;
    }

    @Override
    public void modelRestored(Model model) {
        substitutionProcessKnown = false;
    }

    @Override
    public String toString() {
        if (COUNT_TOTAL_OPERATIONS) {
            return "\tsimulateCount = " + simulateCount + "\n" +
                    "\tgetTraitCount = " + getTraitCount + "\n";
        } else {
            return super.toString();
        }
    }

    @Override
    public void simulate(final int[] operations, final int operationCount, final int rootNodeNumber) {
        // We intentionally bypass the parent implementation pattern:
        // no per-node simulation is needed; just ensure caches and compute gradient.
        setupStatistics();

        if (COUNT_TOTAL_OPERATIONS) {
            ++simulateCount;
        }
    }

    protected void constructGradientTrait(Helper treeTraitHelper) {
        treeTraitHelper.addTrait(new TreeTrait.DA() {
            @Override
            public String getTraitName() {
                return getGradientTraitName();
            }

            @Override
            public TreeTrait.Intent getIntent() {
                return Intent.WHOLE_TREE;
            }

            @Override
            public double[] getTrait(Tree tree, NodeRef node) {
                return getGradient(node);
            }
        });
    }

    protected void constructHessianTrait(Helper treeTraitHelper) {
        treeTraitHelper.addTrait(new TreeTrait.DA() {
            @Override
            public String getTraitName() {
                return getHessianTraitName();
            }

            @Override
            public TreeTrait.Intent getIntent() {
                return Intent.WHOLE_TREE;
            }

            @Override
            public double[] getTrait(Tree tree, NodeRef node) {
                return getHessian(tree, node);
            }
        });
    }
}