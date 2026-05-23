package dr.evomodel.treedatalikelihood.discrete.beastBasedDiscreteTreeLikelihood;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.treedatalikelihood.discrete.beastBasedDiscreteTreeLikelihood.representations.PreOrderRepresentation;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;

import java.util.ArrayDeque;
import java.util.Arrays;
import java.util.Objects;
import java.util.Queue;
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
    private final double[] tmpSiblingScales;
    private final double[] effectiveBranchLengths;
    private final boolean cacheOnlyBranchTopPreOrder;
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

        this.tmpSiblingScales = new double[patternCount];

        this.preOrderAtBranchStart = cacheBranchStartPreOrder ? new double[nodeCount][nodeBufferLength] : null;
        this.preOrderAtBranchEnd = cacheBranchEndPreOrder ? new double[nodeCount][nodeBufferLength] : null;

        this.nodePreOrderKnown = new boolean[nodeCount];
        this.preOrderStartKnown = preOrderAtBranchStart != null ? new boolean[nodeCount] : null;
        this.preOrderEndKnown = preOrderAtBranchEnd != null ? new boolean[nodeCount] : null;

        this.tmpParentNodePreOrder = new double[stateCount];

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

        final Queue<NodeRef> queue = new ArrayDeque<>();
        queue.add(root);

        while (!queue.isEmpty()) {
            final NodeRef parent = queue.remove();
            if (tree.isExternal(parent)) {
                continue;
            }

            final NodeRef left = tree.getChild(parent, 0);
            final NodeRef right = tree.getChild(parent, 1);

            propagateToChild(parent, left, right, categoryRates);
            propagateToChild(parent, right, left, categoryRates);

            queue.add(left);
            queue.add(right);
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

    private void propagateToChild(NodeRef parent,
                                  NodeRef child,
                                  NodeRef sibling,
                                  double[] categoryRates) {
        final int parentNumber = parent.getNumber();
        final int childNumber = child.getNumber();
        if (nodePreOrderKnown[childNumber]) {
            return;
        }

        final int siblingNumber = sibling.getNumber();

        final double[] childPreOrderStart = branchStartPreOrder[childNumber];
        final double[] childPreOrderEnd = cacheOnlyBranchTopPreOrder ? null : nodePreOrder[childNumber];
        final double[] childPreOrderEndStandard =
                nodePreOrderStandard == null
                        ? null
                        : nodePreOrderStandard[childNumber];

        final double[] childScale = nodePreOrderLogScales[childNumber];
        final double[] siblingScale = tmpSiblingScales;
        postOrderMessageProvider.getPostOrderBranchScalesInto(siblingNumber, siblingScale);

        final double childLength = effectiveBranchLengths[childNumber];

        Arrays.fill(childPreOrderStart, 0.0);
        if (childPreOrderEnd != null) {
            Arrays.fill(childPreOrderEnd, 0.0);
        }
        if (childPreOrderEndStandard != null) {
            Arrays.fill(childPreOrderEndStandard, 0.0);
        }

        for (int c = 0; c < categoryCount; c++) {
//            final double siblingEffectiveLength = siblingLength * categoryRates[c];
            final double childEffectiveLength = childLength * categoryRates[c];

            for (int p = 0; p < patternCount; p++) {
                final int off = offset(c, p, 0);
                final double parentNodeScale = prepareParentNodePreOrder(
                        parent, c, p, categoryRates, tmpParentNodePreOrder);
                final double[] parentPreOrder = cacheOnlyBranchTopPreOrder
                        ? (tree.isRoot(parent) ? branchStartPreOrder[parentNumber] : tmpParentNodePreOrder)
                        : (tree.isRoot(parent) ? branchStartPreOrder[parentNumber] : nodePreOrder[parentNumber]);
                final int parentOff = cacheOnlyBranchTopPreOrder && !tree.isRoot(parent) ? 0 : off;

                if (nodePreOrderStandard == null) {
                    final double[] siblingBranchTopPostOrder =
                            postOrderMessageProvider.getPostOrderBranchTopBuffer(siblingNumber);
                    preOrderRepresentation.combineParentAndSibling(
                            parentPreOrder, parentOff,
                            siblingBranchTopPostOrder, off,
                            childPreOrderStart, off
                    );
                } else {
                    final double[] parentNodePreOrderStandard = nodePreOrderStandard[parentNumber];
                    final double[] siblingBranchTopPostOrderStandard =
                            postOrderMessageProvider.getPostOrderBranchTopStandardBuffer(siblingNumber);
                    preOrderRepresentation.importPreOrderProductFromStandard(
                            parentNodePreOrderStandard, parentOff,
                            siblingBranchTopPostOrderStandard, off,
                            childPreOrderStart, off);
                }

                childScale[p] = parentNodeScale + siblingScale[p];

                final double extraScaleStart = normalizePatternSlice(childPreOrderStart, off);
                childScale[p] += extraScaleStart;

                if (!cacheOnlyBranchTopPreOrder) {
                    // propagate down the child branch — writes directly to childPreOrderEnd at off
                    preOrderRepresentation.propagateToBranchBottom(
                            childNumber,
                            childEffectiveLength,
                            childPreOrderStart, off,
                            childPreOrderEnd, off
                    );

                    if (childPreOrderEndStandard != null) {
                        childScale[p] += preOrderRepresentation.normalizeAndExportPreOrderPartialToStandard(
                                childPreOrderEnd, off,
                                childPreOrderEndStandard, off,
                                DEFAULT_SCALING_FLOOR, DEFAULT_SCALING_CEILING);
                    } else {
                        childScale[p] += normalizePatternSlice(childPreOrderEnd, off);
                    }
                }
            }
        }

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
                                             int category,
                                             int pattern,
                                             double[] categoryRates,
                                             double[] outParentNodePreOrder) {
        final int parentNumber = parent.getNumber();
        final int off = offset(category, pattern, 0);
        double scale = nodePreOrderLogScales[parentNumber][pattern];

        if (tree.isRoot(parent)) {
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
