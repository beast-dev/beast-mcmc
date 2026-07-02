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

public abstract class AbstractContinuousBranchRateGradient implements GradientWrtParameterProvider, Reportable {

    protected static final double ZERO_BRANCH_LENGTH_TOLERANCE = 1.0e-12;

    protected final TreeDataLikelihood treeDataLikelihood;
    protected final ContinuousDataLikelihoodDelegate likelihoodDelegate;
    protected final TreeTrait<List<BranchSufficientStatistics>> treeTraitProvider;
    protected final Tree tree;
    protected final Parameter rateParameter;
    protected final DifferentiableBranchRates branchRateModel;
    protected final int nTraits;
    protected Double numericGradientStepSize = null;

    private final String reportName;

    protected AbstractContinuousBranchRateGradient(String reportName,
                                                   String traitName,
                                                   TreeDataLikelihood treeDataLikelihood,
                                                   ContinuousDataLikelihoodDelegate likelihoodDelegate,
                                                   Parameter rateParameter) {

        if (treeDataLikelihood == null) {
            throw new IllegalArgumentException("treeDataLikelihood must be non-null");
        }
        if (likelihoodDelegate == null) {
            throw new IllegalArgumentException("likelihoodDelegate must be non-null");
        }
        if (rateParameter == null) {
            throw new IllegalArgumentException("rateParameter must be non-null");
        }

        final BranchRateModel brm = treeDataLikelihood.getBranchRateModel();
        if (!(brm instanceof DifferentiableBranchRates)) {
            throw new IllegalArgumentException(reportName + " requires differentiable branch rates.");
        }

        this.reportName = reportName;
        this.treeDataLikelihood = treeDataLikelihood;
        this.likelihoodDelegate = likelihoodDelegate;
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

            final int destinationIndex = getParameterIndexFromNode(node);
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

            final double timeDerivative = getBranchTimeDerivative(node, differential);
            if (timeDerivative == 0.0) {
                continue;
            }

            for (int trait = 0; trait < nTraits; ++trait) {
                result[destinationIndex] +=
                        getGradientWrtBranchTime(statisticsForNode.get(trait), node) * timeDerivative;
            }
        }

        return result;
    }

    protected abstract double getGradientWrtBranchTime(BranchSufficientStatistics statistics, NodeRef node);

    protected double getBranchTimeDerivative(NodeRef node, double branchRateDifferential) {
        final double branchRate = branchRateModel.getBranchRate(tree, node);
        if (branchRate == 0.0) {
            return 0.0;
        }
        return getScaledBranchTime(node) * branchRateDifferential / branchRate;
    }

    protected double getScaledBranchTime(NodeRef node) {
        final int matrixIndex = likelihoodDelegate.getDiffusionProcessDelegate()
                .getMatrixBufferOffsetIndex(node.getNumber());
        return likelihoodDelegate.getIntegrator().getBranchLength(matrixIndex);
    }

    protected int getParameterIndexFromNode(NodeRef node) {
        return branchRateModel.getParameterIndexFromNode(node);
    }

    protected boolean isZeroLengthBranch(NodeRef node) {
        return Math.abs(tree.getBranchLength(node)) <= ZERO_BRANCH_LENGTH_TOLERANCE;
    }

    @Override
    public double getNumericGradientStepSize() {
        if (numericGradientStepSize == null) {
            return StepSizeLevel.SMALL.getStepSizeRatio();
        } else {
            return numericGradientStepSize.doubleValue();
        }
    }

    @Override
    public void setNumericGradientStepSize(double stepSize) {
        numericGradientStepSize = stepSize;
    }

    @Override
    public String getReport() {
        return reportName + ".\n" +
                "analytic: " + Arrays.toString(getGradientLogDensity()) + '\n';
    }
}
