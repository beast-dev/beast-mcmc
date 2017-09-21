package dr.evomodel.treedatalikelihood.preorder;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.continuous.MultivariateDiffusionModel;
import dr.evomodel.treedatalikelihood.SimulationTreeTraversal;
import dr.evomodel.treedatalikelihood.continuous.ConjugateRootTraitPrior;
import dr.evomodel.treedatalikelihood.continuous.ContinuousDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.continuous.ContinuousRateTransformation;
import dr.evomodel.treedatalikelihood.continuous.ContinuousTraitPartialsProvider;
import dr.evomodel.treedatalikelihood.continuous.cdi.ContinuousDiffusionIntegrator;
import dr.math.matrixAlgebra.WrappedVector;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;

import java.util.List;

import static dr.math.matrixAlgebra.missingData.MissingOps.unwrap;
import static dr.math.matrixAlgebra.missingData.MissingOps.wrap;

/**
 * @author Marc A. Suchard
 */
public class TipFullConditionalDistributionDelegate extends ProcessSimulationDelegate.AbstractContinuousTraitDelegate {

    public TipFullConditionalDistributionDelegate(String name, Tree tree,
                                                  MultivariateDiffusionModel diffusionModel,
                                                  ContinuousTraitPartialsProvider dataModel,
                                                  ConjugateRootTraitPrior rootPrior,
                                                  ContinuousRateTransformation rateTransformation,
                                                  BranchRateModel rateModel,
                                                  ContinuousDataLikelihoodDelegate likelihoodDelegate) {

        super(name, tree, diffusionModel, dataModel, rootPrior, rateTransformation, rateModel, likelihoodDelegate);
//            buffer = new MeanAndVariance[tree.getExternalNodeCount()];
        this.likelihoodDelegate = likelihoodDelegate;
        this.cdi = likelihoodDelegate.getIntegrator();

        this.dimPartial = dimTrait + likelihoodDelegate.getPrecisionType().getMatrixLength(dimTrait);
        partialNodeBuffer = new double[numTraits * dimPartial];

    }

    protected boolean isLoggable() {
        return false;
    }

    protected final ContinuousDataLikelihoodDelegate likelihoodDelegate;
    protected final ContinuousDiffusionIntegrator cdi;
    protected final int dimPartial;
    protected final double[] partialNodeBuffer;

    public static String getTraitName(String name) {
        return "fcd." + name;
    }

    protected String delegateGetTraitName() {
        return getTraitName(name);
    }

    protected Class delegateGetTraitClass() {
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
                return getTrait(tree, node).toString();
            }

            public boolean getLoggable() {
                return isLoggable();
            }
        };

        treeTraitHelper.addTrait(baseTrait);
    }

    protected double[] getTraitForNode(NodeRef node) {
//        private double[] getTraitForNode(NodeRef node) {

        assert simulationProcess != null;
        assert node != null;

        simulationProcess.cacheSimulatedTraits(node);

        double[] partial = new double[dimPartial * numTraits];
        cdi.getPreOrderPartial(likelihoodDelegate.getActiveNodeIndex(node.getNumber()), partial);

        return partial;
//            return new MeanAndVariance(partial);
    }

    @Override
    protected void simulateRoot(int rootNumber) {

        if (DEBUG) {
            System.err.println("computeRoot");
        }

        final DenseMatrix64F diffusion = new DenseMatrix64F(likelihoodDelegate.getDiffusionModel().getPrecisionmatrix());

        // Copy post-order root prior to pre-order

        final double[] priorBuffer = partialNodeBuffer;
        final double[] rootBuffer = new double[priorBuffer.length];

        cdi.getPostOrderPartial(rootProcessDelegate.getPriorBufferIndex(), partialNodeBuffer); // No double-buffering

        int offset = 0;
        for (int trait = 0; trait < numTraits; ++trait) {
            // Copy mean
            System.arraycopy(priorBuffer, offset, rootBuffer, offset, dimTrait);

            final DenseMatrix64F Pp = wrap(priorBuffer, offset + dimTrait, dimTrait, dimTrait);
            final DenseMatrix64F Pr = new DenseMatrix64F(dimTrait, dimTrait);
            CommonOps.mult(Pp, diffusion, Pr);

            unwrap(Pr, rootBuffer, offset + dimTrait);

            offset += dimPartial;
        }

        cdi.setPreOrderPartial(likelihoodDelegate.getActiveNodeIndex(rootNumber), rootBuffer);

        if (DEBUG) {
            System.err.println("Root: " + new WrappedVector.Raw(rootBuffer, 0, rootBuffer.length));
            System.err.println("");
        }
    }

    @Override
    protected void simulateNode(BranchNodeOperation operation, double branchNormalization) {
        throw new RuntimeException("Not implemented");
    }

    @Override
    protected void simulateNode(NodeOperation operation) {

        cdi.updatePreOrderPartial(
                likelihoodDelegate.getActiveNodeIndex(operation.getNodeNumber()),
                likelihoodDelegate.getActiveNodeIndex(operation.getLeftChild()),
                likelihoodDelegate.getActiveMatrixIndex(operation.getLeftChild()),
                likelihoodDelegate.getActiveNodeIndex(operation.getRightChild()),
                likelihoodDelegate.getActiveMatrixIndex(operation.getRightChild())
        );
    }

    @Override
    public final void simulate(final SimulationTreeTraversal treeTraversal,
                               final int rootNodeNumber) {


        final List<NodeOperation> nodeOperations = treeTraversal.getNodeOperations();
        setupStatistics();

        simulateRoot(rootNodeNumber);

        for (NodeOperation operation : nodeOperations) {
            simulateNode(operation);
        }

        if (DEBUG) {
            System.err.println("END OF PRE-ORDER");
        }
    }

//        final private MeanAndVariance[] buffer;
//
//        private NodeRef nodeForLastCall = null;
//        private MeanAndVariance cachedMeanAndVariance;

    private static final boolean DEBUG = false;

}
