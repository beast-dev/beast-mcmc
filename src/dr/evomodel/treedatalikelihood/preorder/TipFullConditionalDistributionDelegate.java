package dr.evomodel.treedatalikelihood.preorder;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
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
public class TipFullConditionalDistributionDelegate extends ProcessSimulationDelegate.AbstractContinuousTraitDelegate {

    public TipFullConditionalDistributionDelegate(String name, Tree tree,
                                                  MultivariateDiffusionModel diffusionModel,
                                                  ContinuousTraitPartialsProvider dataModel,
                                                  ConjugateRootTraitPrior rootPrior,
                                                  ContinuousRateTransformation rateTransformation,
//                                                  BranchRateModel rateModel,
                                                  ContinuousDataLikelihoodDelegate likelihoodDelegate) {

        super(name, tree, diffusionModel, dataModel, rootPrior, rateTransformation, likelihoodDelegate);
//            buffer = new MeanAndVariance[tree.getExternalNodeCount()];
        this.likelihoodDelegate = likelihoodDelegate;
        this.cdi = likelihoodDelegate.getIntegrator();

        assert (cdi instanceof ContinuousDiffusionIntegrator.Basic);

        this.dimPartial = dimTrait + likelihoodDelegate.getPrecisionType().getMatrixLength(dimTrait);
        partialNodeBuffer = new double[numTraits * dimPartial];
//        partialPriorBuffer = new double[numTraits * dimPartial];
        partialRootBuffer = new double[numTraits * dimPartial];

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
//    private final double[] partialPriorBuffer;
    private final double[] partialRootBuffer;

    public static String getName(String name) {
        return "fcd." + name;
    }

    public String getTraitName(String name) {
        return getName(name);
    }

    private String delegateGetTraitName() {
        return getTraitName(name);
    }

    private Class delegateGetTraitClass() {
        return double[].class;
    }

    protected void constructTraits(Helper treeTraitHelper) {

        TreeTrait.DA baseTrait = new TreeTrait.DA() {

            public String getTraitName() {
                return delegateGetTraitName();
            }

            public Intent getIntent() {
                return Intent.NODE;
            }

            public Class getTraitClass() {
                return delegateGetTraitClass();
            }

            public double[] getTrait(Tree t, NodeRef node) {
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

    private static String formatted(double[] values) {

        if (values.length == 1) {
            return Double.toString(values[0]);
        }

        StringBuilder sb = new StringBuilder();
        sb.append("{");

        for (int i = 0; i < values.length; ++i) {
            sb.append(Double.toString(values[i]));
            if (i < (values.length - 1)) {
                sb.append(",");
            }
        }

        sb.append("}");
        return sb.toString();
    }

    protected double[] getTraitForNode(NodeRef node) {

        assert simulationProcess != null;
        assert node != null;

        simulationProcess.cacheSimulatedTraits(node);

        double[] partial = new double[dimPartial * numTraits];
        cdi.getPreOrderPartial(likelihoodDelegate.getActiveNodeIndex(node.getNumber()), partial);

        return partial;
    }

    @Override
    protected void simulateRoot(final int rootIndex) {

        if (DEBUG) {
            System.err.println("Simulate root node " + rootIndex);
        }

        // Copy from prior to root pre-order buffer
        cdi.getPostOrderPartial(rootProcessDelegate.getPriorBufferIndex(), partialRootBuffer); // No double-buffering
        cdi.setPreOrderPartial(likelihoodDelegate.getActiveNodeIndex(rootIndex), partialRootBuffer); // With double-buffering

        if (DEBUG) {
            System.err.println("Root: " + new WrappedVector.Raw(partialRootBuffer, 0, partialRootBuffer.length));
            System.err.println("");
        }
    }

//    @Override
//    protected void simulateNode(BranchNodeOperation operation, double branchNormalization) {
//        throw new RuntimeException("Not implemented");
//    }

//    @Override
//    protected void simulateNode(NodeOperation operation) {
//
//        cdi.updatePreOrderPartial(
//                likelihoodDelegate.getActiveNodeIndex(operation.getNodeNumber()),
//                likelihoodDelegate.getActiveNodeIndex(operation.getLeftChild()),
//                likelihoodDelegate.getActiveMatrixIndex(operation.getLeftChild()),
//                likelihoodDelegate.getActiveNodeIndex(operation.getRightChild()),
//                likelihoodDelegate.getActiveMatrixIndex(operation.getRightChild())
//        );
//
//        if (DEBUG) {
//            cdi.getPreOrderPartial(likelihoodDelegate.getActiveNodeIndex(operation.getLeftChild()), partialRootBuffer);
//            System.err.println("Node: "
//                    + operation.getLeftChild() + " "
//                    + new WrappedVector.Raw(partialRootBuffer, 0, partialRootBuffer.length));
//            System.err.println("");
//        }
//    }

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

//    @Override
//    public final void simulate(final SimulationTreeTraversal treeTraversal,
//                               final int rootNodeNumber) {
//
//        final List<NodeOperation> nodeOperations = treeTraversal.getNodeOperations();
//        setupStatistics();
//
//        simulateRoot(rootNodeNumber);
//
//        for (NodeOperation operation : nodeOperations) {
//            simulateNode(operation);
//        }
//
//        if (DEBUG) {
//            System.err.println("END OF PRE-ORDER");
//        }
//    }

    private static final boolean DEBUG = false;

}
