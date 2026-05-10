package dr.evomodel.treedatalikelihood.discrete.beastBasedDiscreteTreeLikelihood.representations;

import dr.evomodel.branchmodel.RewardsAwareBranchModel;

import java.util.Arrays;

/***
 * Conventions:
 *   post-order : p_top    = P(t)   p_bottom
 *   pre-order  : q_bottom = P(t)^T q_top
 *
 * In this representation, "internal" and "standard" coordinates coincide,
 * so most conversion methods are just copies.
 */
/*
 * @author Filippo Monti
 */
public final class RewardsAwarePartialsRepresentation
        implements BidirectionalRepresentation {

    private final RewardsAwareBranchModel branchModel;
    private final int stateCount;

    public RewardsAwarePartialsRepresentation(
            RewardsAwareBranchModel branchModel
    ) {
        if (branchModel == null) {
            throw new IllegalArgumentException("branchModel must be non-null");
        }
        this.branchModel = branchModel;
        this.stateCount = branchModel.getStateCount();
    }

    @Override
    public String getName() {
        return "standard";
    }

    @Override
    public int getStateCount() {
        return stateCount;
    }

    @Override
    public void markDirty() {
    }

    @Override
    public void storeState() {
    }

    @Override
    public void restoreState() {
    }

    @Override
    public void updateForLikelihood() {
    }

    @Override
    public void initializeRootPartial(double[] rootFrequencies, double[] outRootPreOrder) {
        checkLength(rootFrequencies, "rootFrequencies");
        checkLength(outRootPreOrder, "outRootPreOrder");

        System.arraycopy(rootFrequencies, 0, outRootPreOrder, 0, stateCount);
    }

    @Override
    public void combineParentAndSibling(
            double[] parentNodePreOrder,
            double[] siblingBranchTopPostOrderStandard,
            double[] outChildBranchTopPreOrder
    ) {
        checkLength(parentNodePreOrder, "parentNodePreOrder");
        checkLength(siblingBranchTopPostOrderStandard, "siblingBranchTopPostOrderStandard");
        checkLength(outChildBranchTopPreOrder, "outChildBranchTopPreOrder");

        for (int i = 0; i < stateCount; i++) {
            outChildBranchTopPreOrder[i] =
                    parentNodePreOrder[i] * siblingBranchTopPostOrderStandard[i];
        }
        boolean allZero = true;
        for (int i = 0; i < stateCount; i++) {
            if (outChildBranchTopPreOrder[i] != 0.0) {
                allZero = false;
                break;
            }
        }

//        if (allZero) {
//            throw new IllegalStateException(
//                    "combineParentAndSibling produced all-zero preorder.\n" +
//                            "parentNodePreOrder=" + Arrays.toString(parentNodePreOrder) + "\n" +
//                            "siblingBranchTopPostOrderStandard=" + Arrays.toString(siblingBranchTopPostOrderStandard)
//            );
//        }
    }

    @Override
    public void propagateToBranchBottom(
            int childNodeNumber,
            double branchLength,
            double[] childBranchTopPreOrder,
            double[] outChildNodePreOrder
    ) {
        checkLength(childBranchTopPreOrder, "childBranchTopPreOrder");
        checkLength(outChildNodePreOrder, "outChildNodePreOrder");

        applyTranspose(
                childNodeNumber,
                childBranchTopPreOrder,
                outChildNodePreOrder
        );
    }

    @Override
    public void initializeTipPartial(double[] standardTip, double[] outPartial) {
        checkLength(standardTip, "standardTip");
        checkLength(outPartial, "outPartial");

        System.arraycopy(standardTip, 0, outPartial, 0, stateCount);
    }

    @Override
    public void propagateToBranchTop(
            int nodeNumber,
            double branchLength,
            double[] childPartial,
            double[] outBranchTopPartial
    ) {
        checkLength(childPartial, "childPartial");
        checkLength(outBranchTopPartial, "outBranchTopPartial");

        apply(nodeNumber,
                childPartial,
                outBranchTopPartial
        );
    }

    @Override
    public void combineBranchTopPartials(
            double[] leftBranchTopPartial,
            double[] rightBranchTopPartial,
            double[] outParentPartial
    ) {
        checkLength(leftBranchTopPartial, "leftBranchTopPartial");
        checkLength(rightBranchTopPartial, "rightBranchTopPartial");
        checkLength(outParentPartial, "outParentPartial");

        for (int i = 0; i < stateCount; i++) {
            outParentPartial[i] =
                    leftBranchTopPartial[i] * rightBranchTopPartial[i];
        }
    }

    private void applyAtomic(int nodeNum, double[] input, double[] output) {
        Arrays.fill(output, 0.0);
        final int stateIndex = branchModel.getAtomicBranchState(nodeNum);
        final double scale = branchModel.getAtomicBranchScale(nodeNum);
        output[stateIndex] = scale * input[stateIndex];
//        int nonZeroCount = 0;
//        int nonZeroIndex = -1;
//        for (int i = 0; i < output.length; i++) {
//            if (output[i] != 0.0) {
//                nonZeroCount++;
//                nonZeroIndex = i;
//            }
//        }
//
//        if (nonZeroCount == 0) {
//            throw new IllegalStateException(
//                    "Atomic branch " + nodeNum +
//                            " produced no nonzero entries in preorder output."
//            );
//        }
//
//        if (nonZeroCount > 1) {
//            throw new IllegalStateException(
//                    "Atomic branch " + nodeNum +
//                            " produced more than one nonzero entry in preorder output."
//            );
//        }
//
//        if (nonZeroCount == 1 && nonZeroIndex != stateIndex) {
//            throw new IllegalStateException(
//                    "Atomic branch " + nodeNum +
//                            " produced a nonzero entry at state " + nonZeroIndex +
//                            " instead of expected state " + stateIndex + "."
//            );
//        }
    }

    public void apply(int nodeNum, double[] input, double[] output) {
        // branchLength is intentionally ignored because RewardsAwareBranchModel
        // reads the current branch length directly from the tree when building
        // transitions.

        if (branchModel.isAtomicBranch(nodeNum)) {
            applyAtomic(nodeNum, input, output);
            return;
        }

        final double[] matrix = branchModel.getTransitionMatrix(nodeNum);

        Arrays.fill(output, 0.0);
        int k = 0;
        for (int i = 0; i < stateCount; i++) {
            double sum = 0.0;
            for (int j = 0; j < stateCount; j++) {
                sum += matrix[k++] * input[j];
            }
            output[i] = sum;
        }
    }

    public void applyTranspose(int nodeNum, double[] input, double[] output) {
        // branchLength is intentionally ignored because RewardsAwareBranchModel
        // reads the current branch length directly from the tree when building
        // transitions.

        if (branchModel.isAtomicBranch(nodeNum)) {
            applyAtomic(nodeNum, input, output);
            return;
        }

        final double[] matrix = branchModel.getTransitionMatrix(nodeNum);

        Arrays.fill(output, 0.0);
        for (int i = 0; i < stateCount; i++) {
            final double in = input[i];
            final int rowOffset = i * stateCount;
            for (int j = 0; j < stateCount; j++) {
                output[j] += matrix[rowOffset + j] * in;
            }
        }
        int nonZeroCount = 0;
        for (int i = 0; i < output.length; i++) {
            if (output[i] != 0.0) {
                nonZeroCount++;
            }
        }

        if (nonZeroCount == 0) {
            throw new IllegalStateException(
                    "Atomic branch " + nodeNum +
                            " produced no nonzero entries in preorder output."
            );
        }
    }

    @Override
    public double rootContribution(double[] rootFrequencies, double[] rootPartial) {
        checkLength(rootFrequencies, "rootFrequencies");
        checkLength(rootPartial, "rootPartial");

        double sum = 0.0;
        for (int i = 0; i < stateCount; i++) {
            sum += rootFrequencies[i] * rootPartial[i];
        }
        return sum;
    }

    @Override
    public void exportPostOrderPartial(double[] partial, double[] outPartial) {
        checkLength(partial, "partial");
        checkLength(outPartial, "outPartial");

        System.arraycopy(partial, 0, outPartial, 0, stateCount);
    }

    @Override
    public void exportPreOrderPartial(double[] partial, double[] outPartial) {
        checkLength(partial, "partial");
        checkLength(outPartial, "outPartial");

        System.arraycopy(partial, 0, outPartial, 0, stateCount);
    }

    private void checkLength(double[] x, String name) {
        if (x == null) {
            throw new IllegalArgumentException(name + " must be non-null");
        }
        if (x.length < stateCount) {
            throw new IllegalArgumentException(
                    name + " has length " + x.length + " but expected at least " + stateCount
            );
        }
    }
}
