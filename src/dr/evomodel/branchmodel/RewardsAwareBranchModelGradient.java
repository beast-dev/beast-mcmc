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

    private final double[] prePartial;
    private final double[] postPartial;
    private final double[] differential;
    private final double[] prePartialAtNode;

    private final double[] gradientBuffer;

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

        this.prePartial = new double[nstates];
        this.prePartialAtNode = new double[nstates];
        this.postPartial = new double[nstates];
        this.differential = new double[nstates * nstates];

        this.gradientBuffer = new double[parameter.getDimension()];
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
        treeDataLikelihood.calculatePostOrderStatistics();
//        discreteDelegate.invalidatePreOrderOnlyForDebug();
        discreteDelegate.ensurePreOrderComputed();
//        treeDataLikelihood.getLogLikelihood();

        Arrays.fill(prePartial, 0.0);
        Arrays.fill(prePartialAtNode, 0.0);
        Arrays.fill(postPartial, 0.0);
        Arrays.fill(gradientBuffer, 0.0);

        for (int b = 0; b < parameter.getDimension(); b++) {
            if (indicator != null && indicator.getParameterValue(b) == indicatorOnBase) {
                gradientBuffer[b] = 0.0;
                continue;
            }

            final int nodeNum = totalRewardsBranchRates.getNodeNumberFromParameterIndex(b);
            final NodeRef node = tree.getNode(nodeNum);

            computeDifferentialForBranch(node, differential);

            discreteDelegate.getPreOrderAtBranchStartInto(nodeNum, prePartial);
            discreteDelegate.getPreOrderAtBranchEndInto(nodeNum, prePartialAtNode);
            discreteDelegate.getPostOrderAtBranchEndInto(nodeNum, postPartial);

//            System.out.println("Branch " + b + " (node " + nodeNum + "): prePartial=" + Arrays.toString(prePartial) + ",\n postPartial=" + Arrays.toString(postPartial) + ",\n differential=" + Arrays.toString(differential));

            final double acc = bilinearFormStable(prePartial, differential, postPartial, nstates);

            double L = 0.0;
            for (int i = 0; i < nstates; i++) {
                L += prePartialAtNode[i] * postPartial[i];
            }
            gradientBuffer[b] = acc / L;
        }

        return gradientBuffer;
    }

    private void computeDifferentialForBranch(NodeRef node, double[] differential) {
        final double branchLength = tree.getBranchLength(node);
        final double totalReward = totalRewardsBranchRates.getBranchRate(tree, node);
        sericola.computePdfDerivativeWrtRhoInto(totalReward, branchLength, differential, false);
    }

    static double bilinearFormStable(double[] pre, double[] D, double[] post, int n) {
        double acc = 0.0, cAcc = 0.0;

        for (int i = 0; i < n; i++) {
            double pre_i = pre[i];
            if (pre_i == 0.0) continue;

            int rowBase = i * n;
            double rowDot = 0.0, cRow = 0.0;

            for (int j = 0; j < n; j++) {
                double y = D[rowBase + j] * post[j] - cRow;
                double t = rowDot + y;
                cRow = (t - rowDot) - y;
                rowDot = t;
            }

            double y = pre_i * rowDot - cAcc;
            double t = acc + y;
            cAcc = (t - acc) - y;
            acc = t;
        }
        return acc;
    }

    @Override
    public String getReport() {
        if (tolerance == null) {
            return "Reward-aware branch-model gradient (no check tolerance specified).";
        } else {
            return "Reward-aware branch-model gradient; check tolerance=" + tolerance;
        }
    }

    public static String getGradientTraitName(String traitName) {
        return GRADIENT_TRAIT_PREFIX + ":" + traitName;
    }

    public static String getHessianTraitName(String traitName) {
        return HESSIAN_TRAIT_PREFIX + ":" + traitName;
    }
}
