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
    private final double[] tmpParentNodePreOrderStandard;
    private final double[] tmpSiblingPostOrder;
    private final double[] tmpSiblingPostOrderStandard;
    private final double[] tmpChildBranchTopPreOrderStandard;
    private final double[] tmpChildBranchTopPreOrder;
    private final double[] tmpChildNodePreOrder;
    private final double[] tmpExport;

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
        this.nodePreOrder = new double[nodeCount][nodeBufferLength];
        this.storedNodePreOrder = new double[nodeCount][nodeBufferLength];
        if (preOrderRepresentation.storesPartialsInStandardBasis()) {
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
        this.tmpParentNodePreOrderStandard = new double[stateCount];
        this.tmpSiblingPostOrder = new double[stateCount];
        this.tmpSiblingPostOrderStandard = new double[stateCount];
        this.tmpChildBranchTopPreOrderStandard = new double[stateCount];
        this.tmpChildBranchTopPreOrder = new double[stateCount];
        this.tmpChildNodePreOrder = new double[stateCount];
        this.tmpExport = new double[stateCount];

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
        copy2D(nodePreOrder, storedNodePreOrder);
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

        double[][] tmpA = nodePreOrder;
        nodePreOrder = storedNodePreOrder;
        storedNodePreOrder = tmpA;

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
        final double[] rootBuffer = nodePreOrder[rootNumber];
        final double[] rootStartBuffer = branchStartPreOrder[rootNumber];
        final double[] rootStandardBuffer = nodePreOrderStandard == null ? null : nodePreOrderStandard[rootNumber];
        Arrays.fill(rootBuffer, 0.0);
        Arrays.fill(rootStartBuffer, 0.0);
        if (rootStandardBuffer != null) {
            Arrays.fill(rootStandardBuffer, 0.0);
        }
        Arrays.fill(nodePreOrderLogScales[rootNumber], 0.0);

        for (int c = 0; c < categoryCount; c++) {
            for (int p = 0; p < patternCount; p++) {
                final int off = offset(c, p, 0);
                preOrderRepresentation.initializeRootPartial(rootFrequencies, tmpChildNodePreOrder);
                System.arraycopy(tmpChildNodePreOrder, 0, rootBuffer, off, stateCount);
                if (rootStandardBuffer != null) {
                    System.arraycopy(rootFrequencies, 0, rootStandardBuffer, off, stateCount);
                }
            }
        }
        System.arraycopy(rootBuffer, 0, rootStartBuffer, 0, rootBuffer.length);

        nodePreOrderKnown[rootNumber] = true;

        // Root has a branch-start preorder value, but exported caches are filled lazily.
    }

//    private void initializeRoot(int rootNumber, double[] rootFrequencies) {
//        final double[] rootBuffer = nodePreOrder[rootNumber];
//        Arrays.fill(rootBuffer, 0.0);
//        Arrays.fill(nodePreOrderLogScales[rootNumber], 0.0);
//
//        for (int c = 0; c < categoryCount; c++) {
//            for (int p = 0; p < patternCount; p++) {
//                final int off = offset(c, p, 0);
//                preOrderRepresentation.initializeRootPartial(rootFrequencies, tmpChildNodePreOrder);
//                System.arraycopy(tmpChildNodePreOrder, 0, rootBuffer, off, stateCount);
//            }
//        }
//
//        nodePreOrderKnown[rootNumber] = true;
//
//        if (preOrderAtBranchEnd != null) {
//            exportNode(rootNumber, preOrderAtBranchEnd[rootNumber]);
//            preOrderEndKnown[rootNumber] = true;
//        }
//    }

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

        // Root uses root-start preorder.
        // Non-root parent uses the END of the parent branch, i.e. preorder at the node.
        final double[] parentNodePreOrder =
                tree.isRoot(parent)
                        ? branchStartPreOrder[parentNumber]
                        : nodePreOrder[parentNumber];
        final double[] parentNodePreOrderStandard =
                nodePreOrderStandard == null
                        ? null
                        : nodePreOrderStandard[parentNumber];

        final double[] childPreOrderStart = branchStartPreOrder[childNumber];
        final double[] childPreOrderEnd = nodePreOrder[childNumber];
        final double[] childPreOrderEndStandard =
                nodePreOrderStandard == null
                        ? null
                        : nodePreOrderStandard[childNumber];

        final double[] childScale = nodePreOrderLogScales[childNumber];
        final double[] parentScale = nodePreOrderLogScales[parentNumber];
        final double[] siblingScale = tmpSiblingScales;
        postOrderMessageProvider.getPostOrderBranchScalesInto(siblingNumber, siblingScale);

        final double childLength = effectiveBranchLengths[childNumber];

        Arrays.fill(childPreOrderStart, 0.0);
        Arrays.fill(childPreOrderEnd, 0.0);
        if (childPreOrderEndStandard != null) {
            Arrays.fill(childPreOrderEndStandard, 0.0);
        }

        for (int c = 0; c < categoryCount; c++) {
//            final double siblingEffectiveLength = siblingLength * categoryRates[c];
            final double childEffectiveLength = childLength * categoryRates[c];

            for (int p = 0; p < patternCount; p++) {
                final int off = offset(c, p, 0);

                if (nodePreOrderStandard == null) {
                    // sibling branch-top postorder in the paired post-order representation
                    postOrderMessageProvider.getPostOrderBranchTopInto(siblingNumber, c, p, tmpSiblingPostOrder);

                    // parent preorder at the correct side, in the representation's internal basis
                    System.arraycopy(parentNodePreOrder, off, tmpParentNodePreOrder, 0, stateCount);

                    // combine parent preorder with sibling branch-top postorder
                    preOrderRepresentation.combineParentAndSibling(
                            tmpParentNodePreOrder,
                            tmpSiblingPostOrder,
                            tmpChildBranchTopPreOrder
                    );
                } else {
                    postOrderMessageProvider.getPostOrderBranchTopStandardInto(
                            siblingNumber, c, p, tmpSiblingPostOrderStandard);
                    System.arraycopy(parentNodePreOrderStandard, off,
                            tmpParentNodePreOrderStandard, 0, stateCount);
                    for (int s = 0; s < stateCount; s++) {
                        tmpChildBranchTopPreOrderStandard[s] =
                                tmpParentNodePreOrderStandard[s] * tmpSiblingPostOrderStandard[s];
                    }
                    preOrderRepresentation.importPreOrderPartialFromStandard(
                            tmpChildBranchTopPreOrderStandard, tmpChildBranchTopPreOrder);
                }

                System.arraycopy(tmpChildBranchTopPreOrder, 0, childPreOrderStart, off, stateCount);

                // THIS is missing in your new code
                childScale[p] = parentScale[p] + siblingScale[p];

                // normalize START before transpose propagation
                final double extraScaleStart = normalizePatternSlice(childPreOrderStart, off);
                childScale[p] += extraScaleStart;

                // use normalized START vector for transpose propagation
                System.arraycopy(childPreOrderStart, off, tmpChildBranchTopPreOrder, 0, stateCount);

                // propagate down the child branch
                preOrderRepresentation.propagateToBranchBottom(
                        childNumber,
                        childEffectiveLength,
                        tmpChildBranchTopPreOrder,
                        tmpChildNodePreOrder
                );

                System.arraycopy(tmpChildNodePreOrder, 0, childPreOrderEnd, off, stateCount);

                // normalize END too
                final double extraScaleEnd = normalizePatternSlice(childPreOrderEnd, off);
                childScale[p] += extraScaleEnd;

                if (childPreOrderEndStandard != null) {
                    System.arraycopy(childPreOrderEnd, off, tmpChildNodePreOrder, 0, stateCount);
                    preOrderRepresentation.exportPreOrderPartialToStandard(
                            tmpChildNodePreOrder, tmpChildBranchTopPreOrderStandard);
                    System.arraycopy(tmpChildBranchTopPreOrderStandard, 0,
                            childPreOrderEndStandard, off, stateCount);
                }
            }
        }

        nodePreOrderKnown[childNumber] = true;

        for (int i = 0; i < childPreOrderEnd.length; i++) {
            if (!Double.isFinite(childPreOrderEnd[i])) {
                throw new IllegalStateException("Non-finite preorder partial at node " + childNumber
                        + ": " + Arrays.toString(childPreOrderEnd));
            }
        }

    }

    private static final double DEFAULT_SCALING_FLOOR = 1.0e-200;
    private static final double DEFAULT_SCALING_CEILING = 1.0e200;

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

    public double[] getPreOrderAtBranchStart(int nodeNumber) {
        final double[] out = new double[flattenedLength()];
        getPreOrderAtBranchStartInto(nodeNumber, out);
        return out;
    }

    public double[] getPreOrderAtBranchEnd(int nodeNumber) {
        final double[] out = new double[flattenedLength()];
        getPreOrderAtBranchEndInto(nodeNumber, out);
        return out;
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
        } else {
            exportInternalBuffer(nodePreOrder[nodeNumber], out);
        }
    }

    public void getPreOrderAtBranchStartInto(int nodeNumber, int category, int pattern, double[] out) {
        final int off = offset(category, pattern, 0);
        if (preOrderAtBranchStart != null) {
            ensurePreOrderAtBranchStartExported(nodeNumber);
            System.arraycopy(preOrderAtBranchStart[nodeNumber], off, out, 0, stateCount);
        } else {
            exportInternalSlice(branchStartPreOrder[nodeNumber], off, out);
        }
    }

    public void getPreOrderAtBranchEndInto(int nodeNumber, int category, int pattern, double[] out) {
        final int off = offset(category, pattern, 0);
        if (preOrderAtBranchEnd != null) {
            ensurePreOrderAtBranchEndExported(nodeNumber);
            System.arraycopy(preOrderAtBranchEnd[nodeNumber], off, out, 0, stateCount);
        } else if (tree.isRoot(tree.getNode(nodeNumber))) {
            Arrays.fill(out, 0.0);
        } else {
            exportInternalSlice(nodePreOrder[nodeNumber], off, out);
        }
    }

    public void getInternalPreOrderAtBranchStartInto(int nodeNumber, int category, int pattern, double[] out) {
        final int off = offset(category, pattern, 0);
        System.arraycopy(branchStartPreOrder[nodeNumber], off, out, 0, stateCount);
    }

    public void getInternalPreOrderAtNodeInto(int nodeNumber, int category, int pattern, double[] out) {
        final int off = offset(category, pattern, 0);
        System.arraycopy(nodePreOrder[nodeNumber], off, out, 0, stateCount);
    }

    public double[] getPreOrderBranchScales(int nodeNumber) {
        return Arrays.copyOf(nodePreOrderLogScales[nodeNumber], nodePreOrderLogScales[nodeNumber].length);
    }

//    public double[] getPreOrderBranchScalesInto(int nodeNumber, double[] out) {
//        return Arrays.copyOf(nodePreOrderLogScales[nodeNumber], nodePreOrderLogScales[nodeNumber].length);
//    }
    public void getPreOrderBranchScalesInto(int nodeNumber, double[] out) {
        System.arraycopy(nodePreOrderLogScales[nodeNumber], 0, out, 0, patternCount);
    }

    private void exportInternalBuffer(double[] source, double[] dest) {
        if (dest.length < flattenedLength()) {
            throw new IllegalArgumentException("Destination length must be at least " + flattenedLength());
        }
        for (int c = 0; c < categoryCount; c++) {
            for (int p = 0; p < patternCount; p++) {
                final int off = offset(c, p, 0);
                exportInternalSlice(source, off, tmpExport);
                System.arraycopy(tmpExport, 0, dest, off, stateCount);
            }
        }
    }

    private void exportInternalSlice(double[] source, int off, double[] dest) {
        sliceInto(source, off, tmpParentNodePreOrder);
        preOrderRepresentation.exportPreOrderPartial(tmpParentNodePreOrder, dest);
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

    private void sliceInto(double[] array, int off, double[] dest) {
        System.arraycopy(array, off, dest, 0, stateCount);
    }

    private static void copy2D(double[][] src, double[][] dst) {
        for (int i = 0; i < src.length; i++) {
            System.arraycopy(src[i], 0, dst[i], 0, src[i].length);
        }
    }

    private void copyFullNodeBuffer(double[] source, double[] dest) {
        final int length = flattenedLength();
        if (dest.length < length) {
            throw new IllegalArgumentException("Destination length must be at least " + length);
        }
        System.arraycopy(source, 0, dest, 0, length);
    }

    private static void requireCache(Object cache, String name) {
        if (cache == null) {
            throw new IllegalStateException("Cache not enabled: " + name);
        }
    }
}
