/*
 * MomentBMContinuousTreeLikelihoodEngine.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 *  BEAST is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 *
 */

package dr.evomodel.treedatalikelihood.continuousNew;

import dr.evomodel.treedatalikelihood.ProcessOnTreeDelegate.BranchOperation;
import dr.evomodel.treedatalikelihood.ProcessOnTreeDelegate.NodeOperation;
import dr.evomodel.treedatalikelihood.continuous.ConjugateRootTraitPrior;

import java.util.Arrays;
import java.util.List;

/**
 * Lightweight homogeneous BM pruning backend in the legacy moment form.
 *
 * <p>This class deliberately does not use the full canonical Gaussian state.
 * For homogeneous BM with scalar tip precision, each post-order message is a
 * trait mean vector plus one scalar precision. Branch propagation is therefore
 * scalar variance inflation and node combination is a scalar-weighted average.
 */
final class MomentBMContinuousTreeLikelihoodEngine implements ContinuousTreeLikelihoodEngine {

    private static final double LOG_SQRT_2_PI = 0.5 * Math.log(2.0 * Math.PI);

    private final ContinuousDataLikelihoodDelegateNew delegate;
    private final ConjugateRootTraitPrior rootPrior;
    private final BMProcessModel processModel;
    private final int nodeCount;
    private final int tipCount;
    private final int dim;
    private final int dim2;
    private final int partialDim;
    private final int partialBufferCount;
    private final double halfDim;
    private final double dimLogSqrtTwoPi;
    private final double[] partials;
    private final double[] remainders;
    private final double[] branchTimes;
    private final double[] storedBranchTimes;
    private final int[] storedBranchTimeIndices;
    private final int[] branchTimeStoreEpochs;
    private final double[] traitPrecision;
    private final double[] storedTraitPrecision;
    private final double[] difference;

    private double logDetTraitPrecision;
    private double storedLogDetTraitPrecision;
    private double logPrecisionConstant;
    private double storedLogPrecisionConstant;
    private boolean precisionDirty = true;
    private boolean storedPrecisionDirty = true;
    private int storedBranchTimeCount = 0;
    private int branchTimeStoreEpoch = 1;

    MomentBMContinuousTreeLikelihoodEngine(final ContinuousDataLikelihoodDelegateNew delegate) {
        this.delegate = delegate;
        this.rootPrior = delegate.getRootPrior();
        this.processModel = new BMProcessModel("momentBm.process", delegate.getDiffusionModel());
        this.nodeCount = delegate.getNodeCount();
        this.tipCount = delegate.getTipCount();
        this.dim = delegate.getTraitDim();
        this.dim2 = dim * dim;
        this.partialDim = dim + 1;
        this.partialBufferCount = delegate.getPartialBufferCount();
        this.halfDim = 0.5 * dim;
        this.dimLogSqrtTwoPi = dim * LOG_SQRT_2_PI;
        this.partials = new double[partialBufferCount * partialDim];
        this.remainders = new double[partialBufferCount];
        this.branchTimes = new double[nodeCount];
        this.storedBranchTimes = new double[nodeCount];
        this.storedBranchTimeIndices = new int[nodeCount];
        this.branchTimeStoreEpochs = new int[nodeCount];
        this.traitPrecision = new double[dim2];
        this.storedTraitPrecision = new double[dim2];
        this.difference = new double[dim];

        delegate.loadAllMomentBMTipData(this);
    }

    @Override
    public double calculateLikelihood(final List<BranchOperation> branchOperations,
                                      final List<NodeOperation> nodeOperations,
                                      final int rootNodeNumber) {
        if (delegate.isComputingWishartStatistics()) {
            throw new UnsupportedOperationException(
                    "Moment-form BM likelihood does not yet provide Wishart statistics; use implementation=\"cdi\".");
        }

        if (precisionDirty || delegate.consumeMomentBMDiffusionModelUpdate()) {
            refreshPrecision();
        }

        delegate.syncMomentBMTipData(this);
        delegate.updateMomentBMBranchNormalization();
        for (BranchOperation operation : branchOperations) {
            final int branchNumber = operation.getBranchNumber();
            setBranchTime(branchNumber, delegate.getEffectiveBranchTime(operation));
        }

        delegate.flipMomentBMNodeBuffers(nodeOperations);
        for (NodeOperation operation : nodeOperations) {
            updateNode(operation.getNodeNumber(), operation.getLeftChild(), operation.getRightChild());
        }

        return calculateRootLogLikelihood(rootNodeNumber);
    }

    void setTipPartial(final int tipIndex, final double[] tipPartial) {
        final int tipBuffer = delegate.getActiveNodeIndex(tipIndex);
        final int partialOffset = tipBuffer * partialDim;
        System.arraycopy(tipPartial, 0, partials, partialOffset, dim);
        partials[partialOffset + dim] = tipPartial[dim];
        remainders[tipBuffer] = 0.0;
    }

    @Override
    public void storeState() {
        storedBranchTimeCount = 0;
        branchTimeStoreEpoch++;
        if (branchTimeStoreEpoch == Integer.MAX_VALUE) {
            Arrays.fill(branchTimeStoreEpochs, 0);
            branchTimeStoreEpoch = 1;
        }
        System.arraycopy(traitPrecision, 0, storedTraitPrecision, 0, traitPrecision.length);
        storedLogDetTraitPrecision = logDetTraitPrecision;
        storedLogPrecisionConstant = logPrecisionConstant;
        storedPrecisionDirty = precisionDirty;
    }

    @Override
    public void restoreState() {
        for (int i = 0; i < storedBranchTimeCount; i++) {
            branchTimes[storedBranchTimeIndices[i]] = storedBranchTimes[i];
        }
        storedBranchTimeCount = 0;
        System.arraycopy(storedTraitPrecision, 0, traitPrecision, 0, traitPrecision.length);
        logDetTraitPrecision = storedLogDetTraitPrecision;
        logPrecisionConstant = storedLogPrecisionConstant;
        precisionDirty = storedPrecisionDirty;
    }

    @Override
    public void acceptState() {
        storedBranchTimeCount = 0;
    }

    private void refreshPrecision() {
        processModel.makeDirty();
        processModel.fillTraitPrecisionFlat(traitPrecision);
        logDetTraitPrecision = processModel.getTraitPrecisionLogDeterminant();
        logPrecisionConstant = -dimLogSqrtTwoPi + 0.5 * logDetTraitPrecision;
        precisionDirty = false;
    }

    private void setBranchTime(final int branchNumber, final double branchTime) {
        if (branchTimeStoreEpochs[branchNumber] != branchTimeStoreEpoch) {
            branchTimeStoreEpochs[branchNumber] = branchTimeStoreEpoch;
            storedBranchTimeIndices[storedBranchTimeCount] = branchNumber;
            storedBranchTimes[storedBranchTimeCount] = branchTimes[branchNumber];
            storedBranchTimeCount++;
        }
        branchTimes[branchNumber] = branchTime;
    }

    private void updateNode(final int nodeNumber, final int leftChild, final int rightChild) {
        final int parentBuffer = delegate.getActiveNodeIndex(nodeNumber);
        final int leftBuffer = delegate.getActiveNodeIndex(leftChild);
        final int rightBuffer = delegate.getActiveNodeIndex(rightChild);
        final int parentOffset = parentBuffer * partialDim;
        final int leftOffset = leftBuffer * partialDim;
        final int rightOffset = rightBuffer * partialDim;

        final double leftPrecision = partials[leftOffset + dim];
        final double rightPrecision = partials[rightOffset + dim];
        final double leftBranchPrecision = branchPrecision(leftPrecision, branchTimes[leftChild]);
        final double rightBranchPrecision = branchPrecision(rightPrecision, branchTimes[rightChild]);

        final double parentPrecision = leftBranchPrecision + rightBranchPrecision;
        partials[parentOffset + dim] = parentPrecision;
        if (parentPrecision == 0.0) {
            Arrays.fill(partials, parentOffset, parentOffset + dim, 0.0);
            remainders[parentBuffer] = remainders[leftBuffer] + remainders[rightBuffer];
            return;
        }

        final double leftWeight = leftBranchPrecision / parentPrecision;
        final double rightWeight = rightBranchPrecision / parentPrecision;

        for (int d = 0; d < dim; d++) {
            partials[parentOffset + d] = leftWeight * partials[leftOffset + d]
                    + rightWeight * partials[rightOffset + d];
        }

        double remainder = remainders[leftBuffer] + remainders[rightBuffer];
        if (leftPrecision != 0.0 && rightPrecision != 0.0) {
            final double remainderPrecision =
                    leftBranchPrecision * rightBranchPrecision / parentPrecision;
            final double ss = quadraticDifference(leftOffset, rightOffset);
            remainder += logPrecisionConstant
                    + halfDim * Math.log(remainderPrecision)
                    - 0.5 * remainderPrecision * ss;
        }
        remainders[parentBuffer] = remainder;
    }

    private double branchPrecision(final double nodePrecision, final double branchTime) {
        return nodePrecision == Double.POSITIVE_INFINITY
                ? 1.0 / branchTime
                : nodePrecision / (1.0 + nodePrecision * branchTime);
    }

    private double calculateRootLogLikelihood(final int rootNodeNumber) {
        final double priorPrecision = rootPrior.getPseudoObservations();
        final int rootBuffer = delegate.getActiveNodeIndex(rootNodeNumber);
        final int rootOffset = rootBuffer * partialDim;
        double rootPrecision = partials[rootOffset + dim];
        if (!Double.isInfinite(priorPrecision)) {
            rootPrecision = rootPrecision * priorPrecision / (rootPrecision + priorPrecision);
        }

        final double ss = rootQuadraticDifference(rootOffset);
        return logPrecisionConstant
                + halfDim * Math.log(rootPrecision)
                - 0.5 * rootPrecision * ss
                + remainders[rootBuffer];
    }

    private double quadraticDifference(final int leftOffset, final int rightOffset) {
        for (int d = 0; d < dim; d++) {
            difference[d] = partials[leftOffset + d] - partials[rightOffset + d];
        }
        return quadraticDifference();
    }

    private double rootQuadraticDifference(final int meanOffset) {
        for (int d = 0; d < dim; d++) {
            difference[d] = partials[meanOffset + d] - rootPrior.getMeanParameter().getParameterValue(d);
        }
        return quadraticDifference();
    }

    private double quadraticDifference() {
        double ss = 0.0;
        for (int row = 0; row < dim; row++) {
            final double rowDifference = difference[row];
            final int rowOffset = row * dim;
            ss += traitPrecision[rowOffset + row] * rowDifference * rowDifference;
            for (int col = 0; col < row; col++) {
                ss += (traitPrecision[rowOffset + col] + traitPrecision[col * dim + row])
                        * rowDifference * difference[col];
            }
        }
        return ss;
    }
}
