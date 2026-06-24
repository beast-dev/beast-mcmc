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

package dr.evomodel.treedatalikelihood.continuous;

import dr.evomodel.treedatalikelihood.ProcessOnTreeDelegate.BranchOperation;
import dr.evomodel.treedatalikelihood.ProcessOnTreeDelegate.NodeOperation;

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

    private final ContinuousDataLikelihoodDelegate delegate;
    private final ConjugateRootTraitPrior rootPrior;
    private final BMProcessModel processModel;
    private final int nodeCount;
    private final int tipCount;
    private final int dim;
    private final int dim2;
    private final double halfDim;
    private final double dimLogSqrtTwoPi;
    private final double[] means;
    private final double[] storedMeans;
    private final double[] precisionScalars;
    private final double[] storedPrecisionScalars;
    private final double[] remainders;
    private final double[] storedRemainders;
    private final double[] branchTimes;
    private final double[] storedBranchTimes;
    private final double[] traitPrecision;
    private final double[] difference;

    private double logDetTraitPrecision;
    private double logPrecisionConstant;
    private boolean precisionDirty = true;

    MomentBMContinuousTreeLikelihoodEngine(final ContinuousDataLikelihoodDelegate delegate) {
        this.delegate = delegate;
        this.rootPrior = delegate.getRootPrior();
        this.processModel = new BMProcessModel("momentBm.process", delegate.getDiffusionModel());
        this.nodeCount = delegate.getNodeCount();
        this.tipCount = delegate.getTipCount();
        this.dim = delegate.getTraitDim();
        this.dim2 = dim * dim;
        this.halfDim = 0.5 * dim;
        this.dimLogSqrtTwoPi = dim * LOG_SQRT_2_PI;
        this.means = new double[nodeCount * dim];
        this.storedMeans = new double[nodeCount * dim];
        this.precisionScalars = new double[nodeCount];
        this.storedPrecisionScalars = new double[nodeCount];
        this.remainders = new double[nodeCount];
        this.storedRemainders = new double[nodeCount];
        this.branchTimes = new double[nodeCount];
        this.storedBranchTimes = new double[nodeCount];
        this.traitPrecision = new double[dim2];
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
            branchTimes[branchNumber] = delegate.getEffectiveBranchTime(operation);
        }

        for (NodeOperation operation : nodeOperations) {
            updateNode(operation.getNodeNumber(), operation.getLeftChild(), operation.getRightChild());
        }

        return calculateRootLogLikelihood(rootNodeNumber);
    }

    void setTipPartial(final int tipIndex, final double[] tipPartial) {
        final int meanOffset = tipIndex * dim;
        System.arraycopy(tipPartial, 0, means, meanOffset, dim);
        precisionScalars[tipIndex] = tipPartial[dim];
        remainders[tipIndex] = 0.0;
    }

    @Override
    public void storeState() {
        System.arraycopy(means, 0, storedMeans, 0, means.length);
        System.arraycopy(precisionScalars, 0, storedPrecisionScalars, 0, precisionScalars.length);
        System.arraycopy(remainders, 0, storedRemainders, 0, remainders.length);
        System.arraycopy(branchTimes, 0, storedBranchTimes, 0, branchTimes.length);
    }

    @Override
    public void restoreState() {
        System.arraycopy(storedMeans, 0, means, 0, means.length);
        System.arraycopy(storedPrecisionScalars, 0, precisionScalars, 0, precisionScalars.length);
        System.arraycopy(storedRemainders, 0, remainders, 0, remainders.length);
        System.arraycopy(storedBranchTimes, 0, branchTimes, 0, branchTimes.length);
        precisionDirty = true;
    }

    @Override
    public void acceptState() {
    }

    private void refreshPrecision() {
        processModel.makeDirty();
        processModel.fillTraitPrecisionFlat(traitPrecision);
        logDetTraitPrecision = processModel.getTraitPrecisionLogDeterminant();
        logPrecisionConstant = -dimLogSqrtTwoPi + 0.5 * logDetTraitPrecision;
        precisionDirty = false;
    }

    private void updateNode(final int nodeNumber, final int leftChild, final int rightChild) {
        final int parentOffset = nodeNumber * dim;
        final int leftOffset = leftChild * dim;
        final int rightOffset = rightChild * dim;

        final double leftPrecision = precisionScalars[leftChild];
        final double rightPrecision = precisionScalars[rightChild];
        final double leftBranchPrecision = branchPrecision(leftPrecision, branchTimes[leftChild]);
        final double rightBranchPrecision = branchPrecision(rightPrecision, branchTimes[rightChild]);

        final double parentPrecision = leftBranchPrecision + rightBranchPrecision;
        precisionScalars[nodeNumber] = parentPrecision;
        if (parentPrecision == 0.0) {
            Arrays.fill(means, parentOffset, parentOffset + dim, 0.0);
            remainders[nodeNumber] = remainders[leftChild] + remainders[rightChild];
            return;
        }

        final double leftWeight = leftBranchPrecision / parentPrecision;
        final double rightWeight = rightBranchPrecision / parentPrecision;

        for (int d = 0; d < dim; d++) {
            means[parentOffset + d] = leftWeight * means[leftOffset + d]
                    + rightWeight * means[rightOffset + d];
        }

        double remainder = remainders[leftChild] + remainders[rightChild];
        if (leftPrecision != 0.0 && rightPrecision != 0.0) {
            final double remainderPrecision =
                    leftBranchPrecision * rightBranchPrecision / parentPrecision;
            final double ss = quadraticDifference(leftOffset, rightOffset);
            remainder += logPrecisionConstant
                    + halfDim * Math.log(remainderPrecision)
                    - 0.5 * remainderPrecision * ss;
        }
        remainders[nodeNumber] = remainder;
    }

    private double branchPrecision(final double nodePrecision, final double branchTime) {
        return nodePrecision == Double.POSITIVE_INFINITY
                ? 1.0 / branchTime
                : nodePrecision / (1.0 + nodePrecision * branchTime);
    }

    private double calculateRootLogLikelihood(final int rootNodeNumber) {
        final double priorPrecision = rootPrior.getPseudoObservations();
        double rootPrecision = precisionScalars[rootNodeNumber];
        if (!Double.isInfinite(priorPrecision)) {
            rootPrecision = rootPrecision * priorPrecision / (rootPrecision + priorPrecision);
        }

        final double ss = rootQuadraticDifference(rootNodeNumber * dim);
        return logPrecisionConstant
                + halfDim * Math.log(rootPrecision)
                - 0.5 * rootPrecision * ss
                + remainders[rootNodeNumber];
    }

    private double quadraticDifference(final int leftOffset, final int rightOffset) {
        for (int d = 0; d < dim; d++) {
            difference[d] = means[leftOffset + d] - means[rightOffset + d];
        }
        return quadraticDifference();
    }

    private double rootQuadraticDifference(final int meanOffset) {
        for (int d = 0; d < dim; d++) {
            difference[d] = means[meanOffset + d] - rootPrior.getMeanParameter().getParameterValue(d);
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
