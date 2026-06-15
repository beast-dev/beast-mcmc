package dr.evomodel.branchmodel;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.branchratemodel.ArbitraryBranchRates;
import dr.evomodel.treedatalikelihood.DataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.DiscreteDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.markovjumps.SericolaSeriesMarkovRewardFastModel;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.xml.Reportable;

import java.util.Arrays;

/**
 * Gradient of log tree data likelihood with respect to branch-specific total rewards
 * stored in an ArbitraryBranchRates object.
 *
 * This version uses DiscreteDataLikelihoodDelegate directly and assumes:
 *  - one pattern
 *  - one category
 *
 * @author Filippo Monti
 */
public final class RewardsAwareBranchModelGradient implements GradientWrtParameterProvider, Reportable {

    private final TreeDataLikelihood treeDataLikelihood;
    private final Tree tree;

    private final DiscreteDataLikelihoodDelegate discreteDelegate;
    private final ArbitraryBranchRates totalRewardsBranchRates;
    private final SericolaSeriesMarkovRewardFastModel sericola;
    private final int nstates;

    private final double indicatorOnBase;

    private final Parameter parameter;
    private final Parameter indicator;

    private final Double tolerance;
    private final boolean useHessian;

    private static final String GRADIENT_TRAIT_PREFIX = "RewardAwareBranchModelGradient";
    private static final String HESSIAN_TRAIT_PREFIX  = "RewardAwareBranchModelHessian";

    private final double[] prePartialAtNode;

    private final double[] gradientBuffer;
    private final double[] batchRewardProportions;
    private final double[] batchBranchLengths;
    private final double[][] batchPrePartials;
    private final double[][] batchPostPartials;
    private final double[] batchDerivativeContractions;
    private final double[] batchLikelihoodDenominators;
    private final int[] batchParameterIndices;

    public RewardsAwareBranchModelGradient(TreeDataLikelihood treeDataLikelihood,
                                              RewardsAwareBranchModel rewardsAwareBranchModel,
                                              ArbitraryBranchRates totalRewardsBranchRates,
                                              Parameter indicator,
                                              Double tolerance,
                                              boolean useHessian,
                                              boolean isOneBasedIndicator) {

        this.treeDataLikelihood = treeDataLikelihood;
        this.tree = rewardsAwareBranchModel.getTree();
        this.totalRewardsBranchRates = totalRewardsBranchRates;
        this.sericola = rewardsAwareBranchModel.getSericolaModel();
        this.nstates = rewardsAwareBranchModel.getStateCount();
        this.indicatorOnBase = isOneBasedIndicator ? 1.0 : 0.0;

        this.tolerance = tolerance;
        this.useHessian = useHessian;

        this.parameter = totalRewardsBranchRates.getRateParameter();
        this.indicator = indicator;

        DataLikelihoodDelegate delegate = treeDataLikelihood.getDataLikelihoodDelegate();
        if (!(delegate instanceof DiscreteDataLikelihoodDelegate)) {
            throw new IllegalArgumentException(
                    "RewardsAwareBranchModelGradient requires TreeDataLikelihood to use DiscreteDataLikelihoodDelegate"
            );
        }
        this.discreteDelegate = (DiscreteDataLikelihoodDelegate) delegate;

        this.prePartialAtNode = new double[nstates];

        this.gradientBuffer = new double[parameter.getDimension()];
        this.batchRewardProportions = new double[parameter.getDimension()];
        this.batchBranchLengths = new double[parameter.getDimension()];
        this.batchPrePartials = new double[parameter.getDimension()][];
        this.batchPostPartials = new double[parameter.getDimension()][];
        this.batchDerivativeContractions = new double[parameter.getDimension()];
        this.batchLikelihoodDenominators = new double[parameter.getDimension()];
        this.batchParameterIndices = new int[parameter.getDimension()];
        for (int i = 0; i < parameter.getDimension(); i++) {
            batchPrePartials[i] = new double[nstates];
            batchPostPartials[i] = new double[nstates];
        }
    }

    @Override
    public Likelihood getLikelihood() {
        return treeDataLikelihood;
    }

    @Override
    public Parameter getParameter() {
        return parameter;
    }

    @Override
    public int getDimension() {
        return parameter.getDimension();
    }

    @Override
    public double[] getGradientLogDensity() {

        // Ensure delegate caches are up to date.
//        treeDataLikelihood.calculatePostOrderStatistics();
        discreteDelegate.updatePostOrdersFromTreeDataLikelihood(treeDataLikelihood);
//        discreteDelegate.invalidatePreOrderOnlyForDebug();
        discreteDelegate.ensurePreOrderComputed();
//        treeDataLikelihood.getLogLikelihood();

        Arrays.fill(prePartialAtNode, 0.0);
        Arrays.fill(gradientBuffer, 0.0);

        final int activeCount = collectActiveContinuousBranchContexts();
        if (activeCount == 0) {
            return gradientBuffer;
        }

        sericola.contractPdfDerivativeWrtRewardProportionInto(
                batchRewardProportions,
                batchBranchLengths,
                activeCount,
                batchPrePartials,
                batchPostPartials,
                batchDerivativeContractions,
                false);

        for (int i = 0; i < activeCount; i++) {
            final int b = batchParameterIndices[i];
            gradientBuffer[b] = batchDerivativeContractions[i] / batchLikelihoodDenominators[i];
        }

        return gradientBuffer;
    }

    private int collectActiveContinuousBranchContexts() {
        int activeCount = 0;
        for (int b = 0; b < parameter.getDimension(); b++) {
            if (indicator != null && indicator.getParameterValue(b) == indicatorOnBase) {
                continue;
            }

            final int nodeNum = totalRewardsBranchRates.getNodeNumberFromParameterIndex(b);
            final NodeRef node = tree.getNode(nodeNum);

            batchParameterIndices[activeCount] = b;
            batchBranchLengths[activeCount] = tree.getBranchLength(node);
            batchRewardProportions[activeCount] = totalRewardsBranchRates.getBranchRate(tree, node);

            final double[] prePartial = batchPrePartials[activeCount];
            final double[] postPartial = batchPostPartials[activeCount];
            discreteDelegate.getPreOrderAtBranchStartInto(nodeNum, prePartial);
            discreteDelegate.getPreOrderAtBranchEndInto(nodeNum, prePartialAtNode);
            discreteDelegate.getPostOrderAtBranchEndInto(nodeNum, postPartial);

            double L = 0.0;
            for (int j = 0; j < nstates; j++) {
                L += prePartialAtNode[j] * postPartial[j];
            }
            batchLikelihoodDenominators[activeCount] = L;

            activeCount++;
        }

        return activeCount;
    }

    @Override
    public String getReport() {
        StringBuilder sb = new StringBuilder();
        if (tolerance == null) {
            sb.append("Reward-aware branch-model gradient (no check tolerance specified).");
        } else {
            sb.append("Reward-aware branch-model gradient; check tolerance=").append(tolerance);
        }
        sb.append('\n');
        sb.append("analytic: ").append(Arrays.toString(getGradientLogDensity())).append('\n');
        return sb.toString();
    }

    public static String getGradientTraitName(String traitName) {
        return GRADIENT_TRAIT_PREFIX + ":" + traitName;
    }

    public static String getHessianTraitName(String traitName) {
        return HESSIAN_TRAIT_PREFIX + ":" + traitName;
    }
}
