package dr.evomodel.treedatalikelihood.discrete.beastBasedDiscreteTreeLikelihood.representations;

import dr.evomodel.branchmodel.RewardsAwareBranchModel;

import java.util.Arrays;

/***
 * Conventions:
 *   post-order : p_top    = P(t)   p_bottom
 *   pre-order  : q_bottom = P(t)^T q_top
 *
 * In this representation, "internal" and "standard" coordinates coincide,
 * so most conversion methods are direct identity writes.
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
        initializeRootPartial(rootFrequencies, outRootPreOrder, 0);
    }

    @Override
    public void initializeRootPartial(double[] rootFrequencies, double[] outRootPreOrder, int outOffset) {
        checkLength(rootFrequencies, "rootFrequencies");
        checkLength(outRootPreOrder, outOffset, "outRootPreOrder");

        for (int i = 0; i < stateCount; i++) {
            outRootPreOrder[outOffset + i] = rootFrequencies[i];
        }
    }

    @Override
    public void combineParentAndSibling(
            double[] parentNodePreOrder,
            int parentOffset,
            double[] siblingBranchTopPostOrderStandard,
            int siblingOffset,
            double[] outChildBranchTopPreOrder,
            int outOffset
    ) {
        checkLength(parentNodePreOrder, parentOffset, "parentNodePreOrder");
        checkLength(siblingBranchTopPostOrderStandard, siblingOffset, "siblingBranchTopPostOrderStandard");
        checkLength(outChildBranchTopPreOrder, outOffset, "outChildBranchTopPreOrder");

        for (int i = 0; i < stateCount; i++) {
            outChildBranchTopPreOrder[outOffset + i] =
                    parentNodePreOrder[parentOffset + i] * siblingBranchTopPostOrderStandard[siblingOffset + i];
        }
        boolean allZero = true;
        for (int i = 0; i < stateCount; i++) {
            if (outChildBranchTopPreOrder[outOffset + i] != 0.0) {
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
            int childBranchTopOffset,
            double[] outChildNodePreOrder,
            int outOffset
    ) {
        checkLength(childBranchTopPreOrder, childBranchTopOffset, "childBranchTopPreOrder");
        checkLength(outChildNodePreOrder, outOffset, "outChildNodePreOrder");

        applyTranspose(childNodeNumber, childBranchTopPreOrder, childBranchTopOffset,
                outChildNodePreOrder, outOffset);
    }

    @Override
    public void initializeTipPartial(double[] standardTip, double[] outPartial) {
        initializeTipPartial(standardTip, outPartial, 0);
    }

    @Override
    public void initializeTipPartial(double[] standardTip, double[] outPartial, int outOffset) {
        checkLength(standardTip, "standardTip");
        checkLength(outPartial, outOffset, "outPartial");

        for (int i = 0; i < stateCount; i++) {
            outPartial[outOffset + i] = standardTip[i];
        }
    }

    @Override
    public void propagateToBranchTop(
            int nodeNumber,
            double branchLength,
            double[] childPartial,
            double[] outBranchTopPartial
    ) {
        propagateToBranchTop(nodeNumber, branchLength, childPartial, 0,
                outBranchTopPartial, 0);
    }

    @Override
    public void propagateToBranchTop(
            int nodeNumber,
            double branchLength,
            double[] childPartial,
            int childOffset,
            double[] outBranchTopPartial,
            int outOffset
    ) {
        checkLength(childPartial, childOffset, "childPartial");
        checkLength(outBranchTopPartial, outOffset, "outBranchTopPartial");

        apply(nodeNumber, childPartial, childOffset, outBranchTopPartial, outOffset);
    }

    @Override
    public void combineBranchTopPartials(
            double[] leftBranchTopPartial,
            double[] rightBranchTopPartial,
            double[] outParentPartial
    ) {
        combineBranchTopPartials(leftBranchTopPartial, 0, rightBranchTopPartial, 0,
                outParentPartial, 0);
    }

    @Override
    public void combineBranchTopPartials(
            double[] leftBranchTopPartial,
            int leftOffset,
            double[] rightBranchTopPartial,
            int rightOffset,
            double[] outParentPartial,
            int outOffset
    ) {
        checkLength(leftBranchTopPartial, leftOffset, "leftBranchTopPartial");
        checkLength(rightBranchTopPartial, rightOffset, "rightBranchTopPartial");
        checkLength(outParentPartial, outOffset, "outParentPartial");

        for (int i = 0; i < stateCount; i++) {
            outParentPartial[outOffset + i] =
                    leftBranchTopPartial[leftOffset + i] * rightBranchTopPartial[rightOffset + i];
        }
    }

    private void applyAtomic(int nodeNum, double[] input, int inputOffset, double[] output, int outputOffset) {
        Arrays.fill(output, outputOffset, outputOffset + stateCount, 0.0);
        final int stateIndex = branchModel.getAtomicBranchState(nodeNum);
        final double scale = branchModel.getAtomicBranchScale(nodeNum);
        output[outputOffset + stateIndex] = scale * input[inputOffset + stateIndex];
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
        apply(nodeNum, input, 0, output, 0);
    }

    public void apply(int nodeNum, double[] input, int inputOffset, double[] output, int outputOffset) {
        // branchLength is intentionally ignored because RewardsAwareBranchModel
        // reads the current branch length directly from the tree when building
        // transitions.

        if (branchModel.isAtomicBranch(nodeNum)) {
            applyAtomic(nodeNum, input, inputOffset, output, outputOffset);
            return;
        }

        final double[] matrix = branchModel.getTransitionMatrix(nodeNum);

        Arrays.fill(output, outputOffset, outputOffset + stateCount, 0.0);
        int k = 0;
        for (int i = 0; i < stateCount; i++) {
            double sum = 0.0;
            for (int j = 0; j < stateCount; j++) {
                sum += matrix[k++] * input[inputOffset + j];
            }
            output[outputOffset + i] = sum;
        }
    }

    public void applyTranspose(int nodeNum, double[] input, int inputOffset, double[] output, int outputOffset) {
        // branchLength is intentionally ignored because RewardsAwareBranchModel
        // reads the current branch length directly from the tree when building
        // transitions.

        if (branchModel.isAtomicBranch(nodeNum)) {
            applyAtomic(nodeNum, input, inputOffset, output, outputOffset);
            return;
        }

        final double[] matrix = branchModel.getTransitionMatrix(nodeNum);

        Arrays.fill(output, outputOffset, outputOffset + stateCount, 0.0);
        for (int i = 0; i < stateCount; i++) {
            final double in = input[inputOffset + i];
            final int rowOffset = i * stateCount;
            for (int j = 0; j < stateCount; j++) {
                output[outputOffset + j] += matrix[rowOffset + j] * in;
            }
        }
        int nonZeroCount = 0;
        for (int i = 0; i < stateCount; i++) {
            if (output[outputOffset + i] != 0.0) {
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
        return rootContribution(rootFrequencies, rootPartial, 0);
    }

    @Override
    public double rootContribution(double[] rootFrequencies, double[] rootPartial, int rootOffset) {
        checkLength(rootFrequencies, "rootFrequencies");
        checkLength(rootPartial, rootOffset, "rootPartial");

        double sum = 0.0;
        for (int i = 0; i < stateCount; i++) {
            sum += rootFrequencies[i] * rootPartial[rootOffset + i];
        }
        return sum;
    }

    @Override
    public void exportPostOrderPartial(double[] partial, double[] outPartial) {
        exportPostOrderPartial(partial, 0, outPartial, 0);
    }

    @Override
    public void exportPostOrderPartial(double[] partial, int partialOffset,
                                       double[] outPartial, int outOffset) {
        checkLength(partial, partialOffset, "partial");
        checkLength(outPartial, outOffset, "outPartial");

        for (int i = 0; i < stateCount; i++) {
            outPartial[outOffset + i] = partial[partialOffset + i];
        }
    }

    @Override
    public void exportPreOrderPartial(double[] partial, double[] outPartial) {
        exportPreOrderPartial(partial, 0, outPartial, 0);
    }

    @Override
    public void exportPreOrderPartial(double[] partial, int partialOffset,
                                      double[] outPartial, int outOffset) {
        checkLength(partial, partialOffset, "partial");
        checkLength(outPartial, outOffset, "outPartial");

        for (int i = 0; i < stateCount; i++) {
            outPartial[outOffset + i] = partial[partialOffset + i];
        }
    }

    private void checkLength(double[] x, String name) {
        checkLength(x, 0, name);
    }

    private void checkLength(double[] x, int offset, String name) {
        if (x == null) {
            throw new IllegalArgumentException(name + " must be non-null");
        }
        if (offset < 0 || x.length < offset + stateCount) {
            throw new IllegalArgumentException(
                    name + " has length " + x.length + " but expected at least " + (offset + stateCount)
            );
        }
    }
}
