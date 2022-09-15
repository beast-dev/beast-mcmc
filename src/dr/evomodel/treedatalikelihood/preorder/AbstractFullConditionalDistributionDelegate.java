package dr.evomodel.treedatalikelihood.preorder;

import dr.evolution.tree.Tree;
import dr.evomodel.continuous.MultivariateDiffusionModel;
import dr.evomodel.treedatalikelihood.continuous.ConjugateRootTraitPrior;
import dr.evomodel.treedatalikelihood.continuous.ContinuousDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.continuous.ContinuousRateTransformation;
import dr.evomodel.treedatalikelihood.continuous.ContinuousTraitPartialsProvider;
import dr.evomodel.treedatalikelihood.continuous.cdi.ContinuousDiffusionIntegrator;
import dr.math.matrixAlgebra.WrappedVector;

import java.util.List;

/**
 * @author Marc A. Suchard
 */
public abstract class AbstractFullConditionalDistributionDelegate
        extends ProcessSimulationDelegate.AbstractContinuousTraitDelegate {

    AbstractFullConditionalDistributionDelegate(String name, Tree tree,
                                                MultivariateDiffusionModel diffusionModel,
                                                ContinuousTraitPartialsProvider dataModel,
                                                ConjugateRootTraitPrior rootPrior,
                                                ContinuousRateTransformation rateTransformation,
                                                ContinuousDataLikelihoodDelegate likelihoodDelegate) {

        super(name, tree, diffusionModel, dataModel, rootPrior, rateTransformation, likelihoodDelegate);
        this.likelihoodDelegate = likelihoodDelegate;
        this.cdi = likelihoodDelegate.getIntegrator();
        this.dimPartial = likelihoodDelegate.getPrecisionType().getPartialsDimension(dimTrait);
        this.partialNodeBuffer = new double[numTraits * dimPartial];
        this.partialRootBuffer = new double[numTraits * dimPartial];
    }

    public int vectorizeNodeOperations(final List<NodeOperation> nodeOperations,
                                       final int[] operations) {
        return likelihoodDelegate.vectorizeNodeOperations(nodeOperations, operations);
    }

    protected boolean isLoggable() {
        return false;
    }

    protected final ContinuousDataLikelihoodDelegate likelihoodDelegate;
    protected final ContinuousDiffusionIntegrator cdi;
    protected final int dimPartial;
    final double[] partialNodeBuffer;
    private final double[] partialRootBuffer;

    protected abstract void constructTraits(Helper treeTraitHelper);

    @Override
    protected void simulateRoot(final int rootIndex) {

        if (DEBUG) {
            System.err.println("Simulate root node " + rootIndex);
        }

        cdi.calculatePreOrderRoot(rootProcessDelegate.getPriorBufferIndex(),
                likelihoodDelegate.getActiveNodeIndex(rootIndex), likelihoodDelegate.getActivePrecisionIndex(0));

        if (DEBUG) {
            cdi.getPreOrderPartial(likelihoodDelegate.getActiveNodeIndex(rootIndex), partialRootBuffer);
            System.err.println("Root: " + new WrappedVector.Raw(partialRootBuffer, 0, partialRootBuffer.length));
            System.err.println("");
        }
    }

    @Override
    protected void simulateNode(final int parentNumber,
                                final int nodeNumber,
                                final int nodeMatrix,
                                final int siblingNumber,
                                final int siblingMatrix) {

        cdi.updatePreOrderPartial(
                parentNumber,
                nodeNumber,
                nodeMatrix,
                siblingNumber,
                siblingMatrix);

        if (DEBUG) {
            cdi.getPreOrderPartial(nodeNumber, partialRootBuffer);
            System.err.println("Node: " + nodeNumber + " "
                    + new WrappedVector.Raw(partialRootBuffer, 0, partialRootBuffer.length));
            System.err.println("");
        }
    }

    private static final boolean DEBUG = false;

}
