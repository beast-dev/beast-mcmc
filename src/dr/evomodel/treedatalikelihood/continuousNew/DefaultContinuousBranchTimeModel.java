/*
 * DefaultContinuousBranchTimeModel.java
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
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 *
 */

package dr.evomodel.treedatalikelihood.continuousNew;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.treedatalikelihood.ProcessOnTreeDelegate.BranchOperation;
import dr.evomodel.treedatalikelihood.continuous.ContinuousRateTransformation;

/**
 * One-segment branch-time model matching the current continuous-trait behavior.
 */
public final class DefaultContinuousBranchTimeModel implements ContinuousBranchTimeModel {

    private final Tree tree;
    private final BranchRateModel rateModel;
    private final ContinuousRateTransformation rateTransformation;
    private double normalization;

    public DefaultContinuousBranchTimeModel(final ContinuousRateTransformation rateTransformation) {
        this(null, null, rateTransformation);
    }

    public DefaultContinuousBranchTimeModel(final Tree tree,
                                            final BranchRateModel rateModel,
                                            final ContinuousRateTransformation rateTransformation) {
        this.tree = tree;
        this.rateModel = rateModel;
        this.rateTransformation = rateTransformation;
        this.normalization = 1.0;
    }

    @Override
    public double updateNormalization() {
        normalization = rateTransformation == null ? 1.0 : rateTransformation.getNormalization();
        return normalization;
    }

    @Override
    public double getNormalization() {
        return normalization;
    }

    @Override
    public double getEffectiveBranchTime(final BranchOperation branchOperation) {
        return branchOperation.getBranchLength() * normalization;
    }

    @Override
    public double getEffectiveBranchTime(final int childNodeIndex) {
        requireTree();
        final NodeRef node = tree.getNode(childNodeIndex);
        if (tree.isRoot(node)) {
            return 0.0;
        }
        double branchTime = tree.getBranchLength(node);
        if (rateModel != null) {
            branchTime *= rateModel.getBranchRate(tree, node);
        }
        return branchTime * normalization;
    }

    @Override
    public int getSegmentCount(final int childNodeIndex) {
        requireTree();
        return tree.isRoot(tree.getNode(childNodeIndex)) ? 0 : 1;
    }

    @Override
    public int getSegmentProcessIndex(final int childNodeIndex, final int segmentIndex) {
        requireSingleSegment(childNodeIndex, segmentIndex);
        return 0;
    }

    @Override
    public double getSegmentEffectiveTime(final int childNodeIndex, final int segmentIndex) {
        requireSingleSegment(childNodeIndex, segmentIndex);
        return getEffectiveBranchTime(childNodeIndex);
    }

    private void requireSingleSegment(final int childNodeIndex, final int segmentIndex) {
        if (segmentIndex != 0 || getSegmentCount(childNodeIndex) != 1) {
            throw new IllegalArgumentException("Invalid one-segment branch time request for node "
                    + childNodeIndex + " segment " + segmentIndex);
        }
    }

    private void requireTree() {
        if (tree == null) {
            throw new UnsupportedOperationException("Tree-indexed branch time requires a tree-backed model.");
        }
    }
}
