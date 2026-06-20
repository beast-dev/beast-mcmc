package dr.evomodel.treedatalikelihood.continuous;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.branchratemodel.DifferentiableBranchRates;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.preorder.BranchConditionalDistributionDelegate;
import dr.evomodel.treedatalikelihood.preorder.BranchSufficientStatistics;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.xml.Reportable;

import java.util.Arrays;
import java.util.List;

public final class OUBranchRateGradient implements GradientWrtParameterProvider, Reportable {

    private static final double ZERO_BRANCH_LENGTH_TOLERANCE = 1.0e-12;

    private final TreeDataLikelihood treeDataLikelihood;
    private final TreeTrait<List<BranchSufficientStatistics>> treeTraitProvider;
    private final Tree tree;
    private final Parameter rateParameter;
    private final DifferentiableBranchRates branchRateModel;
    private final OUBranchTimeGradient branchTimeGradient;
    private final int nTraits;

    public OUBranchRateGradient(String traitName,
                                TreeDataLikelihood treeDataLikelihood,
                                ContinuousDataLikelihoodDelegate likelihoodDelegate,
                                Parameter rateParameter) {

        if (treeDataLikelihood == null) {
            throw new IllegalArgumentException("treeDataLikelihood must be non-null");
        }

        final BranchRateModel brm = treeDataLikelihood.getBranchRateModel();
        if (!(brm instanceof DifferentiableBranchRates)) {
            throw new IllegalArgumentException("OUBranchRateGradient requires differentiable branch rates.");
        }

        this.treeDataLikelihood = treeDataLikelihood;
        this.tree = treeDataLikelihood.getTree();
        this.rateParameter = rateParameter;
        this.branchRateModel = (DifferentiableBranchRates) brm;

        final String bcdName = BranchConditionalDistributionDelegate.getName(traitName);
        if (treeDataLikelihood.getTreeTrait(bcdName) == null) {
            likelihoodDelegate.addBranchConditionalDensityTrait(traitName);
        }

        @SuppressWarnings("unchecked")
        final TreeTrait<List<BranchSufficientStatistics>> unchecked =
                treeDataLikelihood.getTreeTrait(bcdName);
        this.treeTraitProvider = unchecked;
        if (treeTraitProvider == null) {
            throw new IllegalStateException("Unable to create branch conditional density trait: " + bcdName);
        }

        this.nTraits = likelihoodDelegate.getTraitCount();
        if (nTraits != 1) {
            throw new RuntimeException("Not yet implemented for >1 traits");
        }

        this.branchTimeGradient = new OUBranchTimeGradient(
                likelihoodDelegate.getTraitDim(),
                tree,
                likelihoodDelegate);
    }

    @Override
    public Likelihood getLikelihood() {
        return treeDataLikelihood;
    }

    @Override
    public Parameter getParameter() {
        return rateParameter;
    }

    @Override
    public int getDimension() {
        return rateParameter.getDimension();
    }

    @Override
    public double[] getGradientLogDensity() {

        treeDataLikelihood.makeDirty();

        final double[] result = new double[getDimension()];

        for (int i = 0; i < tree.getNodeCount(); ++i) {
            final NodeRef node = tree.getNode(i);
            if (tree.isRoot(node) || isZeroLengthBranch(node)) {
                continue;
            }

            final int destinationIndex = branchRateModel.getParameterIndexFromNode(node);
            assert (destinationIndex != -1);

            final double differential = branchRateModel.getBranchRateDifferential(tree, node);
            if (differential == 0.0) {
                continue;
            }

            final List<BranchSufficientStatistics> statisticsForNode = treeTraitProvider.getTrait(tree, node);
            if (statisticsForNode.size() != nTraits) {
                throw new IllegalStateException(
                        "Expected " + nTraits + " branch sufficient statistics but found " +
                                statisticsForNode.size());
            }

            final double timeDerivative = getTimeDerivative(node, differential);
            for (int trait = 0; trait < nTraits; ++trait) {
                result[destinationIndex] +=
                        branchTimeGradient.getGradientForBranch(statisticsForNode.get(trait), node)[0] *
                                timeDerivative;
            }
        }

        return result;
    }

    private double getTimeDerivative(NodeRef node, double branchRateDifferential) {
        final double branchRate = branchRateModel.getBranchRate(tree, node);
        if (branchRate == 0.0) {
            return 0.0;
        }
        final int matrixIndex = ((ContinuousDataLikelihoodDelegate) treeDataLikelihood.getDataLikelihoodDelegate())
                .getDiffusionProcessDelegate()
                .getMatrixBufferOffsetIndex(node.getNumber());
        final double scaledTime = ((ContinuousDataLikelihoodDelegate) treeDataLikelihood.getDataLikelihoodDelegate())
                .getIntegrator()
                .getBranchLength(matrixIndex);
        return scaledTime * branchRateDifferential / branchRate;
    }

    private boolean isZeroLengthBranch(NodeRef node) {
        return Math.abs(tree.getBranchLength(node)) <= ZERO_BRANCH_LENGTH_TOLERANCE;
    }

    @Override
    public String getReport() {
        return "OU branch-rate gradient (branch-time derivative).\n" +
                "analytic: " + Arrays.toString(getGradientLogDensity()) + '\n';
    }
}
