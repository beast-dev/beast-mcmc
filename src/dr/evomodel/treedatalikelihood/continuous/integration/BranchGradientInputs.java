package dr.evomodel.treedatalikelihood.continuous.integration;

import dr.evomodel.continuous.ou.orthogonalblockdiagonal.OrthogonalBlockPreparedBranchBasis;
import dr.evomodel.treedatalikelihood.continuous.framework.CanonicalPreparedBranchSnapshot;
import dr.evomodel.treedatalikelihood.continuous.gaussian.CanonicalGaussianState;
import dr.evomodel.treedatalikelihood.continuous.gaussian.message.CanonicalGaussianMessageOps;
import dr.evomodel.treedatalikelihood.continuous.gaussian.message.CanonicalLocalTransitionAdjoints;

import java.util.Arrays;

/**
 * Frozen branch-local adjoints prepared from current post-order and pre-order passes.
 */
public final class BranchGradientInputs {
    private final int dimension;
    private final int capacity;
    private final int nodeSlotCount;
    private final int[] activeChildIndices;
    private final double[] activeBranchLengths;
    private final CanonicalLocalTransitionAdjoints[] activeAdjoints;
    private final OrthogonalBlockPreparedBranchBasis[] activeOrthogonalPreparedBasis;
    private final boolean[] stagedActiveByChild;
    private final double[] stagedBranchLengthsByChild;
    private final CanonicalLocalTransitionAdjoints[] stagedAdjointsByChild;
    private final OrthogonalBlockPreparedBranchBasis[] stagedOrthogonalPreparedBasisByChild;
    private final CanonicalGaussianState rootPreOrder;
    private final CanonicalGaussianState rootPostOrder;
    private int activeBranchCount;
    private double rootDiffusionScale;
    private boolean hasOrthogonalPreparedBasis;

    BranchGradientInputs(final int capacity, final int dimension) {
        this.dimension = dimension;
        this.capacity = capacity;
        this.nodeSlotCount = capacity + 1;
        this.activeChildIndices = new int[capacity];
        this.activeBranchLengths = new double[capacity];
        this.activeAdjoints = new CanonicalLocalTransitionAdjoints[capacity];
        this.activeOrthogonalPreparedBasis = new OrthogonalBlockPreparedBranchBasis[capacity];
        this.stagedActiveByChild = new boolean[nodeSlotCount];
        this.stagedBranchLengthsByChild = new double[nodeSlotCount];
        this.stagedAdjointsByChild = new CanonicalLocalTransitionAdjoints[nodeSlotCount];
        this.stagedOrthogonalPreparedBasisByChild = new OrthogonalBlockPreparedBranchBasis[nodeSlotCount];
        for (int i = 0; i < capacity; ++i) {
            this.activeAdjoints[i] = new CanonicalLocalTransitionAdjoints(dimension);
        }
        for (int i = 0; i < nodeSlotCount; ++i) {
            this.stagedAdjointsByChild[i] = new CanonicalLocalTransitionAdjoints(dimension);
        }
        this.rootPreOrder = new CanonicalGaussianState(dimension);
        this.rootPostOrder = new CanonicalGaussianState(dimension);
        this.activeBranchCount = 0;
        this.rootDiffusionScale = 0.0;
        this.hasOrthogonalPreparedBasis = false;
    }

    public int getDimension() {
        return dimension;
    }

    public int getActiveBranchCount() {
        return activeBranchCount;
    }

    public int getActiveChildIndex(final int activeIndex) {
        checkActiveIndex(activeIndex);
        return activeChildIndices[activeIndex];
    }

    public double getBranchLength(final int activeIndex) {
        checkActiveIndex(activeIndex);
        return activeBranchLengths[activeIndex];
    }

    public CanonicalLocalTransitionAdjoints getLocalAdjoints(final int activeIndex) {
        checkActiveIndex(activeIndex);
        return activeAdjoints[activeIndex];
    }

    public double getRootDiffusionScale() {
        return rootDiffusionScale;
    }

    public CanonicalGaussianState getRootPreOrderState() {
        return rootPreOrder;
    }

    public CanonicalGaussianState getRootPostOrderState() {
        return rootPostOrder;
    }

    void clear() {
        Arrays.fill(stagedActiveByChild, false);
        activeBranchCount = 0;
        rootDiffusionScale = 0.0;
        hasOrthogonalPreparedBasis = false;
    }

    void checkCompatible(final int expectedCapacity,
                         final int expectedDimension) {
        if (capacity < expectedCapacity) {
            throw new IllegalArgumentException(
                    "BranchGradientInputs capacity " + capacity
                            + " is too small for " + expectedCapacity + " branches.");
        }
        if (dimension != expectedDimension) {
            throw new IllegalArgumentException(
                    "BranchGradientInputs dimension mismatch: "
                            + dimension + " vs " + expectedDimension + ".");
        }
    }

    void addBranch(final int childIndex,
                   final double branchLength,
                   final CanonicalLocalTransitionAdjoints source,
                   final OrthogonalBlockPreparedBranchBasis orthogonalPreparedBasis) {
        if (activeBranchCount >= capacity) {
            throw new IllegalStateException("BranchGradientInputs capacity exceeded.");
        }
        if (orthogonalPreparedBasis != null) {
            hasOrthogonalPreparedBasis = true;
            activeOrthogonalPreparedBasis[activeBranchCount] = orthogonalPreparedBasis;
        }
        activeChildIndices[activeBranchCount] = childIndex;
        activeBranchLengths[activeBranchCount] = branchLength;
        copyAdjoints(source, activeAdjoints[activeBranchCount]);
        activeBranchCount++;
    }

    void addBranch(final CanonicalPreparedBranchSnapshot snapshot,
                   final CanonicalLocalTransitionAdjoints source) {
        addBranch(
                snapshot.getChildNodeIndex(),
                snapshot.getEffectiveBranchLength(),
                source,
                snapshot.getOrthogonalPreparedBasis());
    }

    void stageBranch(final CanonicalPreparedBranchSnapshot snapshot,
                     final CanonicalLocalTransitionAdjoints source) {
        final int childIndex = snapshot.getChildNodeIndex();
        checkChildIndex(childIndex);
        stagedBranchLengthsByChild[childIndex] = snapshot.getEffectiveBranchLength();
        copyAdjoints(source, stagedAdjointsByChild[childIndex]);
        stagedOrthogonalPreparedBasisByChild[childIndex] = snapshot.getOrthogonalPreparedBasis();
        stagedActiveByChild[childIndex] = true;
    }

    void clearStagedBranch(final int childIndex) {
        checkChildIndex(childIndex);
        stagedActiveByChild[childIndex] = false;
    }

    void compactStagedBranches(final int rootIndex,
                               final boolean useOrthogonalPreparedBasis) {
        checkChildIndex(rootIndex);
        activeBranchCount = 0;
        for (int childIndex = 0; childIndex < nodeSlotCount; ++childIndex) {
            if (childIndex == rootIndex || !stagedActiveByChild[childIndex]) {
                continue;
            }
            if (activeBranchCount >= capacity) {
                throw new IllegalStateException("BranchGradientInputs capacity exceeded during compaction.");
            }
            activeChildIndices[activeBranchCount] = childIndex;
            activeBranchLengths[activeBranchCount] = stagedBranchLengthsByChild[childIndex];
            activeAdjoints[activeBranchCount] = stagedAdjointsByChild[childIndex];
            if (useOrthogonalPreparedBasis) {
                activeOrthogonalPreparedBasis[activeBranchCount] =
                        stagedOrthogonalPreparedBasisByChild[childIndex];
            }
            activeBranchCount++;
        }
        hasOrthogonalPreparedBasis = useOrthogonalPreparedBasis && activeBranchCount > 0;
    }

    OrthogonalBlockPreparedBranchBasis getOrthogonalPreparedBasis(final int activeIndex) {
        checkActiveIndex(activeIndex);
        if (!hasOrthogonalPreparedBasis) {
            throw new IllegalStateException("No orthogonal prepared branch basis is available.");
        }
        return activeOrthogonalPreparedBasis[activeIndex];
    }

    void setRoot(final double rootDiffusionScale,
                 final CanonicalGaussianState rootPreOrder,
                 final CanonicalGaussianState rootPostOrder) {
        this.rootDiffusionScale = rootDiffusionScale;
        CanonicalGaussianMessageOps.copyState(rootPreOrder, this.rootPreOrder);
        CanonicalGaussianMessageOps.copyState(rootPostOrder, this.rootPostOrder);
    }

    private void checkActiveIndex(final int activeIndex) {
        if (activeIndex < 0 || activeIndex >= activeBranchCount) {
            throw new IndexOutOfBoundsException(
                    "Active branch index " + activeIndex + " is outside [0, "
                            + activeBranchCount + ").");
        }
    }

    private void checkChildIndex(final int childIndex) {
        if (childIndex < 0 || childIndex >= nodeSlotCount) {
            throw new IndexOutOfBoundsException(
                    "Child index " + childIndex + " is outside [0, " + nodeSlotCount + ").");
        }
    }

    private static void copyAdjoints(final CanonicalLocalTransitionAdjoints source,
                                     final CanonicalLocalTransitionAdjoints target) {
        final int dimension = source.getDimension();
        System.arraycopy(source.dLogL_dF, 0, target.dLogL_dF, 0, dimension * dimension);
        System.arraycopy(source.dLogL_df, 0, target.dLogL_df, 0, dimension);
        System.arraycopy(source.dLogL_dOmega, 0, target.dLogL_dOmega, 0, dimension * dimension);
    }
}
