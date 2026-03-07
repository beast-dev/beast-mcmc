package dr.evomodel.branchmodel;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.branchratemodel.ArbitraryBranchRates;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.preorder.RewardsAwarePartialLikelihoodProvider;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.markovjumps.SericolaSeriesMarkovRewardFastModelRho;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.xml.Reportable;
import java.util.Arrays;

/**
 * Gradient of log tree data likelihood wrt branch-specific total rewards
 * stored in an ArbitraryBranchRates object.
 *
 * This mirrors the overall structure of HomogeneousSubstitutionParameterGradient,
 * but:
 *  - the per-branch differential matrix comes from Sericola (inside RewardsAwareBranchModel)
 *  - the parameter is the rates parameter inside ArbitraryBranchRates
 */

/*
 * @author Filippo Monti
 */
public final class RewardsAwareBranchModelGradient implements GradientWrtParameterProvider, Reportable {

    private final TreeDataLikelihood treeDataLikelihood;
    private final Tree tree;

    private final RewardsAwarePartialLikelihoodProvider partialLikelihoodProvider;
    private final ArbitraryBranchRates totalRewardsBranchRates;
    private final SericolaSeriesMarkovRewardFastModelRho sericola;
    private int nstates;

    private double indicatorOnBase;

    private final Parameter parameter; // the underlying "total.rewards" parameter (dimension = #branches)
    private final Parameter indicator;

    private final Double tolerance;   // optional check tolerance (for Reportable)
    private final boolean useHessian; // reserved; not implemented below

    // ---- Namespacing constants ----
    private static final String GRADIENT_TRAIT_PREFIX = "RewardAwareBranchModelGradient";
    private static final String HESSIAN_TRAIT_PREFIX  = "RewardAwareBranchModelHessian";


    private final double[] prePartial;
    private final double[] postPartial;
    private final double[] differential; // to hold dP/d(totalReward)
    private final double[] prePartialAtNode;

    private final double[] gradientBuffer;


    public RewardsAwareBranchModelGradient(TreeDataLikelihood treeDataLikelihood,
                                           RewardsAwareBranchModel rewardsAwareBranchModel,
                                           RewardsAwarePartialLikelihoodProvider partialLikelihoodProvider,
                                           ArbitraryBranchRates totalRewardsBranchRates,
                                           Parameter indicator,
                                           Double tolerance,
                                           boolean useHessian,
                                           boolean isOneBasedIndicator) {     // if true, indicator=1 means "include in gradient"; if false, indicator=0 means "include in gradient". Only relevant if indicator is provided.

        this.treeDataLikelihood = treeDataLikelihood;
        this.tree = rewardsAwareBranchModel.getTree();
        this.partialLikelihoodProvider = partialLikelihoodProvider;
        this.totalRewardsBranchRates = totalRewardsBranchRates;
        this.sericola = rewardsAwareBranchModel.getSericolaModel();
        this.nstates = rewardsAwareBranchModel.getStateCount();
        this.indicatorOnBase = isOneBasedIndicator ? 1.0 : 0.0;

        this.tolerance = tolerance;
        this.useHessian = useHessian;

        this.parameter = totalRewardsBranchRates.getRateParameter();
        this.indicator = indicator;

        this.prePartial = new double[nstates];
        this.prePartialAtNode = new double[nstates];
        this.postPartial = new double[nstates];
        this.differential = new double[nstates * nstates];

        this.gradientBuffer = new double[ totalRewardsBranchRates.getRateParameter().getDimension() ];
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
        // One dimension per branch-specific total reward
        return parameter.getDimension();
    }

    @Override
    public double[] getGradientLogDensity() {
        partialLikelihoodProvider.updatePreOrderPartials();
        Arrays.fill(prePartial, 0.0);
        Arrays.fill(prePartialAtNode, 0.0);
        Arrays.fill(postPartial, 0.0);
        Arrays.fill(gradientBuffer, 0.0);

        for (int b = 0; b < parameter.getDimension(); b++) {
            if (indicator != null && indicator.getParameterValue(b) == indicatorOnBase) {
                // If indicator is provided, skip branches where indicator=indicatorOnBase
                gradientBuffer[b] = 0;
                continue;
            }
            int nodeNum = totalRewardsBranchRates.getNodeNumberFromParameterIndex(b);
            NodeRef node = tree.getNode(nodeNum);
            computeDifferentialForBranch(node, differential);
            partialLikelihoodProvider.getAbovePartialsForBranch(b, prePartial);
            partialLikelihoodProvider.getAbovePartialsForNode(node, prePartialAtNode);
            partialLikelihoodProvider.getBelowPartialsForBranch(b, postPartial);

            final double acc = bilinearFormStable(prePartial, differential, postPartial, nstates);

            // compute likelihood
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
                // y = D_ij * post_j
                double y = D[rowBase + j] * post[j] - cRow;
                double t = rowDot + y;
                cRow = (t - rowDot) - y;
                rowDot = t;
            }

            // acc += pre_i * rowDot (also Kahan)
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