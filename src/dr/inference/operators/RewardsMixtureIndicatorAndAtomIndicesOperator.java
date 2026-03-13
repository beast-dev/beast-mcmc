package dr.inference.operators;


import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.branchmodel.RewardsAwareBranchModel;
import dr.evomodel.treedatalikelihood.preorder.RewardsAwarePartialLikelihoodProvider;
import dr.inference.model.Parameter;
import dr.math.MathUtils;

import java.util.Arrays;

import static dr.inference.operators.AdaptationMode.ADAPTATION_OFF;

/*
 * RewardsMixtureIndicatorAndAtomIndicesOperator
 *
 * Discrete update for the per-branch mixture indicator z_b ∈ {0,1}
 * and, when z_b = 1, the associated atomic index s_b.
 *
 * For each branch b we compute:
 *
 *   atomic weights:      w1(k) = exp(-λ t_b) * π_k * pre[k] * post[k]
 *   continuous weight:   w0    = preᵀ W_cts post
 *
 * where W_cts is the continuous transition matrix for the branch
 * (already including the continuous-component mass).
 *
 * The regime is sampled using weights
 *
 *   P(z_b = 1) ∝ Σ_k w1(k),    P(z_b = 0) ∝ w0
 *
 * and if z_b = 1 the atomic index is drawn with probabilities
 * proportional to w1(k).
 *
 * pre and post denote the above and below partial likelihood vectors
 * at the branch boundary.
 */

/*
 * @author Filippo Monti
 */

public final class RewardsMixtureIndicatorAndAtomIndicesOperator extends AbstractAdaptableOperator {

    // Inputs
    private final Parameter indicatorZ; // z_b (dimension B), integer-like doubles, values in {0,1}
    private final Parameter atomIndex; // a_b (dimension B), integer-like doubles, values in {0,...,K-1}
    private double updateProportion;
    private final boolean adaptUpdateProportion;
    private final int[] branchBuffer;

    // adaptable value is logit(updateProportion)
    private double adaptableParameter;

    private final RewardsAwareBranchModel rewardsAwareBranchModel;
    private final double[] prePartial;
    private final double[] postPartial;
    private final double[] tmp0;
    private double uniformPrior;
    private final Tree tree;
    private final int nstates;
    private final RewardsAwarePartialLikelihoodProvider partialLikelihoodProvider;

    // No adaptation, but AbstractAdaptableOperator requires an "adaptable parameter".
    private double dummyAdaptable = 1.0;

    public RewardsMixtureIndicatorAndAtomIndicesOperator(
            final Parameter indicatorZ,
            final Parameter atomIndex,
            final Parameter pi,
            final RewardsAwareBranchModel rewardsAwareBranchModel,
            RewardsAwarePartialLikelihoodProvider partialLikelihoodProvider,
            final double updateProportion,
            final boolean adaptUpdateProportion,
            final double weight
    ) {
        super(adaptUpdateProportion ? AdaptationMode.ADAPTATION_ON : ADAPTATION_OFF);


        this.indicatorZ = indicatorZ;
        this.atomIndex = atomIndex;
        this.rewardsAwareBranchModel = rewardsAwareBranchModel;
        this.partialLikelihoodProvider = partialLikelihoodProvider;

        this.tree = rewardsAwareBranchModel.getTree();
        this.updateProportion = updateProportion;
        this.adaptUpdateProportion = adaptUpdateProportion;
        this.adaptableParameter = logit(updateProportion);

        final int B = indicatorZ.getDimension();
        this.branchBuffer = new int[B];
        for (int i = 0; i < B; i++) {
            branchBuffer[i] = i;
        }

        this.nstates = rewardsAwareBranchModel.getStateCount();
        this.prePartial = new double[nstates];
        this.postPartial = new double[nstates];
        this.tmp0 = new double[nstates];

        setWeight(weight);

        if (pi != null && pi.getDimension() != nstates) {
            throw new IllegalArgumentException("pi dimension must equal stateCount K. pi.dim="
                    + pi.getDimension() + ", K=" + nstates);
        } else if (pi == null) {
            // fill piDefault with uniform weights
            this.uniformPrior = 1.0 / nstates;
        }
        if (atomIndex.getDimension() != indicatorZ.getDimension()) {
            throw new IllegalArgumentException("atomIndex and indicatorZ must have the same dimension.");
        }
    }

    @Override
    public String getOperatorName() {
        return "RewardsMixtureIndicatorAndAtomIndicesOperator";
    }

    @Override
    public double doOperation() {

        partialLikelihoodProvider.updatePreOrderPartials();

        final int B = indicatorZ.getDimension();
        final int nToUpdate = Math.max(1, (int) Math.round(updateProportion * B));

        if (nToUpdate >= B) {
            for (int b = 0; b < B; b++) {
                resampleIndicatorForBranch(b);
            }
            return 0.0;
        }

        // partial Fisher-Yates shuffle to sample without replacement
        for (int i = 0; i < B; i++) {
            branchBuffer[i] = i;
        }

        for (int i = 0; i < nToUpdate; i++) {
            final int j = i + MathUtils.nextInt(B - i);
            final int tmp = branchBuffer[i];
            branchBuffer[i] = branchBuffer[j];
            branchBuffer[j] = tmp;

            resampleIndicatorForBranch(branchBuffer[i]);
        }
        atomIndex.fireParameterChangedEvent();
        indicatorZ.fireParameterChangedEvent();

        return 0.0;
    }

    //    @Override
//    public double doOperation() {
//
//        // Ensure pre-order partials are current (child-side buffers as per your convention)
//        partialLikelihoodProvider.updatePreOrderPartials();
//
//        final int B = indicatorZ.getDimension();
//
//        if (updateProportion) {
//            for (int b = 0; b < B; b++) {
//                resampleIndicatorForBranch(b);
//            }
//            return 0.0;
//        }
//
//        final int b = MathUtils.nextInt(B);
//        resampleIndicatorForBranch(b);
//        return 0.0;
//    }

    private void resampleIndicatorForBranch(final int b) {

        final int nodeNum = rewardsAwareBranchModel.getNodeNumberForBranchIndex(b);
        final NodeRef node = tree.getNode(nodeNum);
        if (tree.isRoot(node)) return;

        Arrays.fill(prePartial, 0.0);
        Arrays.fill(postPartial, 0.0);
        partialLikelihoodProvider.getAbovePartialsForBranch(b, prePartial);
        partialLikelihoodProvider.getBelowPartialsForBranch(b, postPartial);

        Arrays.fill(tmp0, 0.0);
        double[] atomicWeightsForBranch = tmp0; // reuse tmp0 for atomic weights
        computeAtomicWeightsForBranchInto(b, prePartial, postPartial, atomicWeightsForBranch);
        double ctsWeight = computeCtsWeightForBranch(b, prePartial, postPartial);
        if (Double.isNaN(ctsWeight) || Double.isInfinite(ctsWeight) || ctsWeight < 0.0) {
            throw new IllegalStateException("Invalid continuous weight: " + ctsWeight);
        }
        double atomicWeight = 0.0;
        for (int j = 0; j < nstates; j++) {
            if (Double.isNaN(atomicWeightsForBranch[j]) || Double.isInfinite(atomicWeightsForBranch[j]) || atomicWeightsForBranch[j] < 0.0) {
                throw new IllegalStateException("Invalid atomic weight at state " + j + ": " + atomicWeightsForBranch[j]);
            }
            atomicWeight += atomicWeightsForBranch[j];
        }

        final double denom = atomicWeight + ctsWeight;
        if (!(denom > 0.0) || Double.isInfinite(denom)) {
            throw new IllegalStateException("Invalid total weight: atomic=" + atomicWeight + ", continuous=" + ctsWeight);
        }
        int newZ = (MathUtils.nextDouble() < (atomicWeight / denom)) ? 1 : 0;
        indicatorZ.setParameterValue(b, newZ);
        if (newZ == 1) resampleAtom(b, atomicWeightsForBranch);
    }

    private void computeAtomicWeightsForBranchInto(final int b,
                                                   double[] prePartial, double[] postPartial,
                                                   double[] atomicWeightsOut) {
        final int nodeNum = rewardsAwareBranchModel.getNodeNumberForBranchIndex(b);
        final NodeRef node = tree.getNode(nodeNum);
        final double uniformizationRate = rewardsAwareBranchModel.getUniformizationRate();
        final double branchLength = tree.getBranchLength(node);
        final double scale = Math.exp(-uniformizationRate * branchLength);
        for (int j = 0; j < nstates; j++) {
            atomicWeightsOut[j] = scale * prePartial[j] * postPartial[j];
        }
    }

    private double computeCtsWeightForBranch(final int b, double[] prePartial, double[] postPartial) {
        final int nodeNum = rewardsAwareBranchModel.getNodeNumberForBranchIndex(b);
        final NodeRef node = tree.getNode(nodeNum);
        if (tree.isRoot(node)) return 0.0;

        final double[] W = rewardsAwareBranchModel.getTransitionMatrixCts(nodeNum);
        return bilinearFormStable(prePartial, W, postPartial, nstates);
    }

    private void resampleAtom(final int b, final double[] atomicWeights) {
        int newIndex = multinomialSampling(atomicWeights);
        atomIndex.setParameterValue(b, newIndex);
    }

    private int multinomialSampling(final double[] atomicWeights) {

        final int nstates = atomicWeights.length;

        // compute total weight
        double total = 0.0;
        for (double w : atomicWeights) {
            total += w;
        }

        if (total <= 0.0) {
            throw new IllegalArgumentException("All atomic weights are zero or negative.");
        }

        // draw uniform in [0, total)
        double u = dr.math.MathUtils.nextDouble() * total;

        // cumulative sampling
        double cumulative = 0.0;
        for (int k = 0; k < nstates; k++) {
            cumulative += atomicWeights[k];
            if (u < cumulative) {
                return k;
            }
        }

        // numerical safety (in case of rounding)
        return nstates - 1;
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

    // ---------------------------------------------------------------------
    // Adaptation
    // ---------------------------------------------------------------------

    @Override
    protected void setAdaptableParameterValue(final double value) {
        if (!adaptUpdateProportion) return;
        adaptableParameter = value;
        updateProportion = logistic(value);

        // keep away from exact 0 and 1
        final double eps = 1e-6;
        if (updateProportion < eps) updateProportion = eps;
        if (updateProportion > 1.0) updateProportion = 1.0;
    }

    @Override
    protected double getAdaptableParameterValue() {
        return adaptableParameter;
    }

    @Override
    public String getAdaptableParameterName() {
        return "updateProportion";
    }

    @Override
    public double getRawParameter() {
        return updateProportion;
    }

    private static double logistic(final double x) {
        if (x >= 0.0) {
            final double e = Math.exp(-x);
            return 1.0 / (1.0 + e);
        } else {
            final double e = Math.exp(x);
            return e / (1.0 + e);
        }
    }

    private static double logit(final double p) {
        final double eps = 1e-12;
        final double pp = Math.max(eps, Math.min(1.0 - eps, p));
        return Math.log(pp / (1.0 - pp));
    }

}