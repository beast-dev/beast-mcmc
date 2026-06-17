package dr.evomodel.treedatalikelihood.discrete.discretetreedataLikelihood;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.treedatalikelihood.discrete.discretetreedataLikelihood.representations.PreOrderRepresentation;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;

import java.util.ArrayDeque;
import java.util.Arrays;
import java.util.Objects;
import java.util.logging.Logger;

/**
 * Pre-order traversal/caching engine for discrete-state models.
 *
 * This class is intentionally separate from the post-order delegate.
 * It depends on post-order only through {@link PostOrderMessageProvider},
 * which exposes sibling messages in the paired post-order representation.
 *
 * Internal pre-order storage is determined by {@link PreOrderRepresentation}.
 */
/*
* @author Filippo Monti
*/
public final class DiscretePreOrderDelegate extends AbstractModel {

    private static final Logger LOGGER = Logger.getLogger("dr.evomodel");

    private final Tree tree;
    private final PreOrderRepresentation preOrderRepresentation;
    private final PostOrderMessageProvider postOrderMessageProvider;
    private final double[] tmpLeftBranchScales;
    private final double[] tmpRightBranchScales;
    private final double[] effectiveBranchLengths;
    private final boolean cacheOnlyBranchTopPreOrder;
    private final ArrayDeque<NodeRef> preOrderQueue;
    private double[] categoryRates;

    private final int nodeCount;
    private final int tipCount;
    private final int patternCount;
    private final int categoryCount;
    private final int stateCount;

    // internal representation buffers: [node][flattened(category,pattern,state)]
    private double[][] branchStartPreOrder;
    private double[][] storedBranchStartPreOrder;
    private double[][] nodePreOrder;
    private double[][] storedNodePreOrder;
    private double[][] nodePreOrderStandard;
    private double[][] storedNodePreOrderStandard;

    // log scales [node][pattern]
    private double[][] nodePreOrderLogScales;
    private double[][] storedNodePreOrderLogScales;

    // optional exported caches
    private final double[][] preOrderAtBranchStart;
    private final double[][] preOrderAtBranchEnd;

    private final boolean[] nodePreOrderKnown;
    private final boolean[] preOrderStartKnown;
    private final boolean[] preOrderEndKnown;

    // scratch
    private final double[] tmpParentNodePreOrder;
    private final double[] tmpPreparedParentPreOrder;

    public DiscretePreOrderDelegate(Tree tree,
                                    double[] branchLengths,
                                    PreOrderRepresentation preOrderRepresentation,
                                    PostOrderMessageProvider postOrderMessageProvider,
                                    int patternCount,
                                    int categoryCount,
                                    boolean cacheBranchStartPreOrder,
                                    boolean cacheBranchEndPreOrder) {
        super("DiscretePreOrderDelegate");

        this.tree = Objects.requireNonNull(tree, "tree");
        this.preOrderRepresentation = Objects.requireNonNull(preOrderRepresentation, "preOrderRepresentation");
        this.postOrderMessageProvider = Objects.requireNonNull(postOrderMessageProvider, "postOrderMessageProvider");
        this.effectiveBranchLengths = Objects.requireNonNull(branchLengths, "branchLengths");
        this.cacheOnlyBranchTopPreOrder = preOrderRepresentation.cacheOnlyBranchTopPreOrder();

        this.nodeCount = tree.getNodeCount();
        this.tipCount = tree.getExternalNodeCount();
        this.patternCount = patternCount;
        this.categoryCount = categoryCount;
        this.stateCount = preOrderRepresentation.getStateCount();

        if (postOrderMessageProvider.getStateCount() != stateCount) {
            throw new IllegalArgumentException("State count mismatch between pre-order and post-order providers");
        }
        if (postOrderMessageProvider.getPatternCount() != patternCount) {
            throw new IllegalArgumentException("Pattern count mismatch between pre-order and post-order providers");
        }
        if (postOrderMessageProvider.getCategoryCount() != categoryCount) {
            throw new IllegalArgumentException("Category count mismatch between pre-order and post-order providers");
        }
        if (categoryCount > 1) {
            throw new IllegalArgumentException("Category count > 1 not currently supported");
        }

        final int nodeBufferLength = flattenedLength();

        this.branchStartPreOrder = new double[nodeCount][nodeBufferLength];
        this.storedBranchStartPreOrder = new double[nodeCount][nodeBufferLength];
        this.nodePreOrder = cacheOnlyBranchTopPreOrder ? null : new double[nodeCount][nodeBufferLength];
        this.storedNodePreOrder = cacheOnlyBranchTopPreOrder ? null : new double[nodeCount][nodeBufferLength];
        if (preOrderRepresentation.storesPartialsInStandardBasis() || cacheOnlyBranchTopPreOrder) {
            this.nodePreOrderStandard = null;
            this.storedNodePreOrderStandard = null;
        } else {
            this.nodePreOrderStandard = new double[nodeCount][nodeBufferLength];
            this.storedNodePreOrderStandard = new double[nodeCount][nodeBufferLength];
        }

        this.nodePreOrderLogScales = new double[nodeCount][patternCount];
        this.storedNodePreOrderLogScales = new double[nodeCount][patternCount];

        this.tmpLeftBranchScales = new double[patternCount];
        this.tmpRightBranchScales = new double[patternCount];
        this.preOrderQueue = new ArrayDeque<NodeRef>(nodeCount);

        this.preOrderAtBranchStart = cacheBranchStartPreOrder ? new double[nodeCount][nodeBufferLength] : null;
        this.preOrderAtBranchEnd = cacheBranchEndPreOrder ? new double[nodeCount][nodeBufferLength] : null;

        this.nodePreOrderKnown = new boolean[nodeCount];
        this.preOrderStartKnown = preOrderAtBranchStart != null ? new boolean[nodeCount] : null;
        this.preOrderEndKnown = preOrderAtBranchEnd != null ? new boolean[nodeCount] : null;

        this.tmpParentNodePreOrder = new double[stateCount];
        this.tmpPreparedParentPreOrder = new double[stateCount];

        LOGGER.info("Creating DiscretePreOrderDelegate");
        LOGGER.info("    Pre-order representation: " + preOrderRepresentation.getName());
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        invalidateAll();
        fireModelChanged();
    }

    @Override
    protected void handleVariableChangedEvent(dr.inference.model.Variable variable, int index,
                                              dr.inference.model.Parameter.ChangeType type) {
        invalidateAll();
        fireModelChanged();
    }

    @Override
    protected void acceptState() {
        // nothing
    }

    public void makeDirty() {
        invalidateAll();
        preOrderRepresentation.markDirty();
        fireModelChanged();
    }

    @Override
    public void storeState() {
        copy2D(branchStartPreOrder, storedBranchStartPreOrder);
        if (nodePreOrder != null) {
            copy2D(nodePreOrder, storedNodePreOrder);
        }
        if (nodePreOrderStandard != null) {
            copy2D(nodePreOrderStandard, storedNodePreOrderStandard);
        }
        copy2D(nodePreOrderLogScales, storedNodePreOrderLogScales);
        preOrderRepresentation.storeState();
    }

    @Override
    public void restoreState() {
        double[][] tmpStart = branchStartPreOrder;
        branchStartPreOrder = storedBranchStartPreOrder;
        storedBranchStartPreOrder = tmpStart;

        if (nodePreOrder != null) {
            double[][] tmpA = nodePreOrder;
            nodePreOrder = storedNodePreOrder;
            storedNodePreOrder = tmpA;
        }

        if (nodePreOrderStandard != null) {
            double[][] tmpStandard = nodePreOrderStandard;
            nodePreOrderStandard = storedNodePreOrderStandard;
            storedNodePreOrderStandard = tmpStandard;
        }

        double[][] tmpB = nodePreOrderLogScales;
        nodePreOrderLogScales = storedNodePreOrderLogScales;
        storedNodePreOrderLogScales = tmpB;

        invalidateAll();
        preOrderRepresentation.restoreState();
    }

    public void ensurePreOrder(int rootNodeNumber,
                               double[] categoryRates,
                               double[] rootFrequencies) {

        this.categoryRates = Objects.requireNonNull(categoryRates, "categoryRates");
        preOrderRepresentation.updateForLikelihood();

        final NodeRef root = tree.getNode(rootNodeNumber);
        final int rootNumber = root.getNumber();

        if (!nodePreOrderKnown[rootNumber]) {
            initializeRoot(rootNumber, rootFrequencies);
        }

        preOrderQueue.clear();
        preOrderQueue.add(root);

        while (!preOrderQueue.isEmpty()) {
            final NodeRef parent = preOrderQueue.remove();
            if (tree.isExternal(parent)) {
                continue;
            }

            final NodeRef left = tree.getChild(parent, 0);
            final NodeRef right = tree.getChild(parent, 1);

            propagateChildren(parent, left, right, categoryRates);

            preOrderQueue.add(left);
            preOrderQueue.add(right);
        }
    }

    private void initializeRoot(int rootNumber, double[] rootFrequencies) {
        final double[] rootBuffer = cacheOnlyBranchTopPreOrder
                ? branchStartPreOrder[rootNumber]
                : nodePreOrder[rootNumber];
        final double[] rootStartBuffer = branchStartPreOrder[rootNumber];
        final double[] rootStandardBuffer = nodePreOrderStandard == null ? null : nodePreOrderStandard[rootNumber];
        Arrays.fill(rootBuffer, 0.0);
        if (rootStartBuffer != rootBuffer) {
            Arrays.fill(rootStartBuffer, 0.0);
        }
        if (rootStandardBuffer != null) {
            Arrays.fill(rootStandardBuffer, 0.0);
        }
        Arrays.fill(nodePreOrderLogScales[rootNumber], 0.0);

        for (int c = 0; c < categoryCount; c++) {
            for (int p = 0; p < patternCount; p++) {
                final int off = offset(c, p, 0);
                preOrderRepresentation.initializeRootPartial(rootFrequencies, rootBuffer, off);
                if (rootStartBuffer != rootBuffer) {
                    preOrderRepresentation.initializeRootPartial(rootFrequencies, rootStartBuffer, off);
                }
                if (rootStandardBuffer != null) {
                    for (int s = 0; s < stateCount; s++) {
                        rootStandardBuffer[off + s] = rootFrequencies[s];
                    }
                }
            }
        }

        nodePreOrderKnown[rootNumber] = true;

        // Root has a branch-start preorder value, but exported caches are filled lazily.
    }

    private void propagateChildren(NodeRef parent,
                                   NodeRef left,
                                   NodeRef right,
                                   double[] categoryRates) {
        final int parentNumber = parent.getNumber();
        final int leftNumber = left.getNumber();
        final int rightNumber = right.getNumber();

        final boolean propagateLeft = !nodePreOrderKnown[leftNumber];
        final boolean propagateRight = !nodePreOrderKnown[rightNumber];
        if (!propagateLeft && !propagateRight) {
            return;
        }

        if (propagateLeft) {
            clearChildPreOrderBuffers(leftNumber);
            postOrderMessageProvider.getPostOrderBranchScalesInto(rightNumber, tmpRightBranchScales);
        }
        if (propagateRight) {
            clearChildPreOrderBuffers(rightNumber);
            postOrderMessageProvider.getPostOrderBranchScalesInto(leftNumber, tmpLeftBranchScales);
        }

        final double[] leftPreOrderStart = branchStartPreOrder[leftNumber];
        final double[] leftPreOrderEnd = cacheOnlyBranchTopPreOrder ? null : nodePreOrder[leftNumber];
        final double[] leftPreOrderEndStandard =
                nodePreOrderStandard == null
                        ? null
                        : nodePreOrderStandard[leftNumber];
        final double[] leftScale = nodePreOrderLogScales[leftNumber];

        final double[] rightPreOrderStart = branchStartPreOrder[rightNumber];
        final double[] rightPreOrderEnd = cacheOnlyBranchTopPreOrder ? null : nodePreOrder[rightNumber];
        final double[] rightPreOrderEndStandard =
                nodePreOrderStandard == null
                        ? null
                        : nodePreOrderStandard[rightNumber];
        final double[] rightScale = nodePreOrderLogScales[rightNumber];

        final double leftLength = effectiveBranchLengths[leftNumber];
        final double rightLength = effectiveBranchLengths[rightNumber];
        final boolean parentIsRoot = tree.isRoot(parent);
        final boolean usePreparedParent =
                propagateLeft && propagateRight
                        && nodePreOrderStandard == null
                        && preOrderRepresentation.supportsPreparedParentForSiblingCombinations();
        for (int c = 0; c < categoryCount; c++) {
            final double leftEffectiveLength = leftLength * categoryRates[c];
            final double rightEffectiveLength = rightLength * categoryRates[c];

            for (int p = 0; p < patternCount; p++) {
                final int off = offset(c, p, 0);
                final double parentNodeScale = prepareParentNodePreOrder(
                        parent, parentIsRoot, c, p, categoryRates, tmpParentNodePreOrder);
                final double[] parentPreOrder = cacheOnlyBranchTopPreOrder
                        ? (parentIsRoot ? branchStartPreOrder[parentNumber] : tmpParentNodePreOrder)
                        : (parentIsRoot ? branchStartPreOrder[parentNumber] : nodePreOrder[parentNumber]);
                final int parentOff = cacheOnlyBranchTopPreOrder && !parentIsRoot ? 0 : off;

                if (usePreparedParent) {
                    preOrderRepresentation.prepareParentForSiblingCombinations(
                            parentPreOrder, parentOff, tmpPreparedParentPreOrder, 0);
                    propagateChildFromPreparedParent(
                            leftNumber, rightNumber,
                            leftPreOrderStart, leftPreOrderEnd, leftPreOrderEndStandard,
                            leftScale, tmpRightBranchScales,
                            leftEffectiveLength, off, p, parentNodeScale);
                    propagateChildFromPreparedParent(
                            rightNumber, leftNumber,
                            rightPreOrderStart, rightPreOrderEnd, rightPreOrderEndStandard,
                            rightScale, tmpLeftBranchScales,
                            rightEffectiveLength, off, p, parentNodeScale);
                } else {
                    if (propagateLeft) {
                        propagateChildFromParent(
                                parentNumber, leftNumber, rightNumber,
                                parentPreOrder, parentOff,
                                leftPreOrderStart, leftPreOrderEnd, leftPreOrderEndStandard,
                                leftScale, tmpRightBranchScales,
                                leftEffectiveLength, off, p, parentNodeScale);
                    }
                    if (propagateRight) {
                        propagateChildFromParent(
                                parentNumber, rightNumber, leftNumber,
                                parentPreOrder, parentOff,
                                rightPreOrderStart, rightPreOrderEnd, rightPreOrderEndStandard,
                                rightScale, tmpLeftBranchScales,
                                rightEffectiveLength, off, p, parentNodeScale);
                    }
                }
            }
        }

        if (propagateLeft) {
            markChildPreOrderKnownAndCheck(leftNumber, leftPreOrderEnd);
        }
        if (propagateRight) {
            markChildPreOrderKnownAndCheck(rightNumber, rightPreOrderEnd);
        }
    }

    private void clearChildPreOrderBuffers(int childNumber) {
        Arrays.fill(branchStartPreOrder[childNumber], 0.0);
        if (!cacheOnlyBranchTopPreOrder) {
            Arrays.fill(nodePreOrder[childNumber], 0.0);
        }
        if (nodePreOrderStandard != null) {
            Arrays.fill(nodePreOrderStandard[childNumber], 0.0);
        }
    }

    private void propagateChildFromPreparedParent(int childNumber,
                                                  int siblingNumber,
                                                  double[] childPreOrderStart,
                                                  double[] childPreOrderEnd,
                                                  double[] childPreOrderEndStandard,
                                                  double[] childScale,
                                                  double[] siblingScale,
                                                  double childEffectiveLength,
                                                  int off,
                                                  int pattern,
                                                  double parentNodeScale) {
        final double[] siblingBranchTopPostOrder =
                postOrderMessageProvider.getPostOrderBranchTopBuffer(siblingNumber);
        preOrderRepresentation.combinePreparedParentAndSibling(
                tmpPreparedParentPreOrder, 0,
                siblingBranchTopPostOrder, off,
                childPreOrderStart, off);
        finishChildPreOrderPattern(
                childNumber,
                childPreOrderStart, childPreOrderEnd, childPreOrderEndStandard,
                childScale, siblingScale,
                childEffectiveLength, off, pattern, parentNodeScale);
    }

    private void propagateChildFromParent(int parentNumber,
                                          int childNumber,
                                          int siblingNumber,
                                          double[] parentPreOrder,
                                          int parentOff,
                                          double[] childPreOrderStart,
                                          double[] childPreOrderEnd,
                                          double[] childPreOrderEndStandard,
                                          double[] childScale,
                                          double[] siblingScale,
                                          double childEffectiveLength,
                                          int off,
                                          int pattern,
                                          double parentNodeScale) {
        if (nodePreOrderStandard == null) {
            final double[] siblingBranchTopPostOrder =
                    postOrderMessageProvider.getPostOrderBranchTopBuffer(siblingNumber);
            preOrderRepresentation.combineParentAndSibling(
                    parentPreOrder, parentOff,
                    siblingBranchTopPostOrder, off,
                    childPreOrderStart, off);
        } else {
            final double[] parentNodePreOrderStandard = nodePreOrderStandard[parentNumber];
            final double[] siblingBranchTopPostOrderStandard =
                    postOrderMessageProvider.getPostOrderBranchTopStandardBuffer(siblingNumber);
            preOrderRepresentation.importPreOrderProductFromStandard(
                    parentNodePreOrderStandard, parentOff,
                    siblingBranchTopPostOrderStandard, off,
                    childPreOrderStart, off);
        }

        finishChildPreOrderPattern(
                childNumber,
                childPreOrderStart, childPreOrderEnd, childPreOrderEndStandard,
                childScale, siblingScale,
                childEffectiveLength, off, pattern, parentNodeScale);
    }

    private void finishChildPreOrderPattern(int childNumber,
                                            double[] childPreOrderStart,
                                            double[] childPreOrderEnd,
                                            double[] childPreOrderEndStandard,
                                            double[] childScale,
                                            double[] siblingScale,
                                            double childEffectiveLength,
                                            int off,
                                            int pattern,
                                            double parentNodeScale) {
        childScale[pattern] = parentNodeScale + siblingScale[pattern];

        final double extraScaleStart = normalizePatternSlice(childPreOrderStart, off);
        childScale[pattern] += extraScaleStart;

        if (!cacheOnlyBranchTopPreOrder) {
            preOrderRepresentation.propagateToBranchBottom(
                    childNumber,
                    childEffectiveLength,
                    childPreOrderStart, off,
                    childPreOrderEnd, off);

            if (childPreOrderEndStandard != null) {
                childScale[pattern] += preOrderRepresentation.normalizeAndExportPreOrderPartialToStandard(
                        childPreOrderEnd, off,
                        childPreOrderEndStandard, off,
                        DEFAULT_SCALING_FLOOR, DEFAULT_SCALING_CEILING);
            } else {
                childScale[pattern] += normalizePatternSlice(childPreOrderEnd, off);
            }
        }
    }

    private void markChildPreOrderKnownAndCheck(int childNumber, double[] childPreOrderEnd) {
        nodePreOrderKnown[childNumber] = true;
        if (childPreOrderEnd != null) {
            for (int i = 0; i < childPreOrderEnd.length; i++) {
                if (!Double.isFinite(childPreOrderEnd[i])) {
                    throw new IllegalStateException("Non-finite preorder partial at node " + childNumber
                            + ": " + Arrays.toString(childPreOrderEnd));
                }
            }
        }

    }

    private static final double DEFAULT_SCALING_FLOOR = 1.0e-200;
    private static final double DEFAULT_SCALING_CEILING = 1.0e200;

    private double prepareParentNodePreOrder(NodeRef parent,
                                             boolean parentIsRoot,
                                             int category,
                                             int pattern,
                                             double[] categoryRates,
                                             double[] outParentNodePreOrder) {
        final int parentNumber = parent.getNumber();
        final int off = offset(category, pattern, 0);
        double scale = nodePreOrderLogScales[parentNumber][pattern];

        if (parentIsRoot) {
            return scale;
        }

        if (!cacheOnlyBranchTopPreOrder) {
            return scale;
        }

        preOrderRepresentation.propagateToBranchBottom(
                parentNumber,
                effectiveBranchLengths[parentNumber] * categoryRates[category],
                branchStartPreOrder[parentNumber], off,
                outParentNodePreOrder, 0);
        scale += normalizePatternSlice(outParentNodePreOrder, 0);
        return scale;
    }

    private double normalizePatternSlice(double[] buffer, int off) {
        double max = 0.0;
        for (int s = 0; s < stateCount; s++) {
            max = Math.max(max, Math.abs(buffer[off + s]));
        }

//        if (max == 0.0) {
//            return Double.NEGATIVE_INFINITY;
//        }
        if (max == 0.0) { // TODO recheck this
            final double uniform = 1.0 / stateCount;
            for (int s = 0; s < stateCount; s++) {
                buffer[off + s] = uniform;
            }
            return Double.NEGATIVE_INFINITY;
        }

        if (max < DEFAULT_SCALING_FLOOR || max > DEFAULT_SCALING_CEILING) {
            for (int s = 0; s < stateCount; s++) {
                buffer[off + s] /= max;
            }
            return Math.log(max);
        }

        return 0.0;
    }

    public void getPreOrderAtBranchStartInto(int nodeNumber, double[] out) {
        if (preOrderAtBranchStart != null) {
            ensurePreOrderAtBranchStartExported(nodeNumber);
            copyFullNodeBuffer(preOrderAtBranchStart[nodeNumber], out);
        } else {
            exportInternalBuffer(branchStartPreOrder[nodeNumber], out);
        }
    }

    public void getPreOrderAtBranchEndInto(int nodeNumber, double[] out) {
        if (preOrderAtBranchEnd != null) {
            ensurePreOrderAtBranchEndExported(nodeNumber);
            copyFullNodeBuffer(preOrderAtBranchEnd[nodeNumber], out);
        } else if (tree.isRoot(tree.getNode(nodeNumber))) {
            Arrays.fill(out, 0, flattenedLength(), 0.0);
        } else if (cacheOnlyBranchTopPreOrder) {
            exportBranchBottomInternalBuffer(nodeNumber, out);
        } else {
            exportInternalBuffer(nodePreOrder[nodeNumber], out);
        }
    }

    public void getPreOrderAtBranchStartInto(int nodeNumber, int category, int pattern, double[] out) {
        final int off = offset(category, pattern, 0);
        if (preOrderAtBranchStart != null) {
            ensurePreOrderAtBranchStartExported(nodeNumber);
            copySliceToOutput(preOrderAtBranchStart[nodeNumber], off, out);
        } else {
            exportInternalSlice(branchStartPreOrder[nodeNumber], off, out);
        }
    }

    public void getPreOrderAtBranchEndInto(int nodeNumber, int category, int pattern, double[] out) {
        final int off = offset(category, pattern, 0);
        if (preOrderAtBranchEnd != null) {
            ensurePreOrderAtBranchEndExported(nodeNumber);
            copySliceToOutput(preOrderAtBranchEnd[nodeNumber], off, out);
        } else if (tree.isRoot(tree.getNode(nodeNumber))) {
            Arrays.fill(out, 0.0);
        } else if (cacheOnlyBranchTopPreOrder) {
            computeInternalPreOrderAtNodeInto(nodeNumber, category, pattern, out, 0);
            preOrderRepresentation.exportPreOrderPartial(out, 0, out, 0);
        } else {
            exportInternalSlice(nodePreOrder[nodeNumber], off, out);
        }
    }

    public void getInternalPreOrderAtBranchStartInto(int nodeNumber, int category, int pattern, double[] out) {
        final int off = offset(category, pattern, 0);
        copySliceToOutput(branchStartPreOrder[nodeNumber], off, out);
    }

    public void getInternalPreOrderAtNodeInto(int nodeNumber, int category, int pattern, double[] out) {
        final int off = offset(category, pattern, 0);
        if (cacheOnlyBranchTopPreOrder) {
            computeInternalPreOrderAtNodeInto(nodeNumber, category, pattern, out, 0);
        } else {
            copySliceToOutput(nodePreOrder[nodeNumber], off, out);
        }
    }

    public void getPreOrderBranchScalesInto(int nodeNumber, double[] out) {
        for (int p = 0; p < patternCount; p++) {
            out[p] = nodePreOrderLogScales[nodeNumber][p];
        }
    }

    private void exportInternalBuffer(double[] source, double[] dest) {
        if (dest.length < flattenedLength()) {
            throw new IllegalArgumentException("Destination length must be at least " + flattenedLength());
        }
        for (int c = 0; c < categoryCount; c++) {
            for (int p = 0; p < patternCount; p++) {
                final int off = offset(c, p, 0);
                preOrderRepresentation.exportPreOrderPartial(source, off, dest, off);
            }
        }
    }

    private void exportBranchBottomInternalBuffer(int nodeNumber, double[] dest) {
        if (dest.length < flattenedLength()) {
            throw new IllegalArgumentException("Destination length must be at least " + flattenedLength());
        }
        for (int c = 0; c < categoryCount; c++) {
            for (int p = 0; p < patternCount; p++) {
                final int off = offset(c, p, 0);
                computeInternalPreOrderAtNodeInto(nodeNumber, c, p, dest, off);
                preOrderRepresentation.exportPreOrderPartial(dest, off, dest, off);
            }
        }
    }

    private void computeInternalPreOrderAtNodeInto(int nodeNumber, int category, int pattern,
                                                   double[] out, int outOffset) {
        if (tree.isRoot(tree.getNode(nodeNumber))) {
            Arrays.fill(out, outOffset, outOffset + stateCount, 0.0);
            return;
        }
        if (categoryRates == null) {
            throw new IllegalStateException("Pre-order has not been computed yet");
        }
        final int off = offset(category, pattern, 0);
        preOrderRepresentation.propagateToBranchBottom(
                nodeNumber,
                effectiveBranchLengths[nodeNumber] * categoryRates[category],
                branchStartPreOrder[nodeNumber], off,
                out, outOffset);
        normalizePatternSlice(out, outOffset);
    }

    private void exportInternalSlice(double[] source, int off, double[] dest) {
        preOrderRepresentation.exportPreOrderPartial(source, off, dest, 0);
    }

    private void ensurePreOrderAtBranchStartExported(int nodeNumber) {
        if (preOrderStartKnown != null && !preOrderStartKnown[nodeNumber]) {
            exportInternalBuffer(branchStartPreOrder[nodeNumber], preOrderAtBranchStart[nodeNumber]);
            preOrderStartKnown[nodeNumber] = true;
        }
    }

    private void ensurePreOrderAtBranchEndExported(int nodeNumber) {
        if (preOrderEndKnown != null && !preOrderEndKnown[nodeNumber]) {
            if (tree.isRoot(tree.getNode(nodeNumber))) {
                Arrays.fill(preOrderAtBranchEnd[nodeNumber], 0.0);
            } else if (cacheOnlyBranchTopPreOrder) {
                exportBranchBottomInternalBuffer(nodeNumber, preOrderAtBranchEnd[nodeNumber]);
            } else {
                exportInternalBuffer(nodePreOrder[nodeNumber], preOrderAtBranchEnd[nodeNumber]);
            }
            preOrderEndKnown[nodeNumber] = true;
        }
    }

    private void invalidateAll() {
        Arrays.fill(nodePreOrderKnown, false);
        if (preOrderStartKnown != null) Arrays.fill(preOrderStartKnown, false);
        if (preOrderEndKnown != null) Arrays.fill(preOrderEndKnown, false);

        for (int i = 0; i < nodeCount; i++) {
            Arrays.fill(nodePreOrderLogScales[i], 0.0);
        }
    }

    private int flattenedLength() {
        return categoryCount * patternCount * stateCount;
    }

    private int offset(int category, int pattern, int state) {
        return ((category * patternCount) + pattern) * stateCount + state;
    }

    private static void copy2D(double[][] src, double[][] dst) {
        for (int i = 0; i < src.length; i++) {
            for (int j = 0; j < src[i].length; j++) {
                dst[i][j] = src[i][j];
            }
        }
    }

    private void copyFullNodeBuffer(double[] source, double[] dest) {
        final int length = flattenedLength();
        if (dest.length < length) {
            throw new IllegalArgumentException("Destination length must be at least " + length);
        }
        for (int i = 0; i < length; i++) {
            dest[i] = source[i];
        }
    }

    private void copySliceToOutput(double[] src, int srcOffset, double[] dst) {
        for (int s = 0; s < stateCount; s++) {
            dst[s] = src[srcOffset + s];
        }
    }
}
