/*
 * CoalescentIntervalTraversal.java
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

package dr.evomodel.coalescent.basta;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.bigfasttree.BigFastTreeIntervals;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.coalescent.StructuredCoalescentActiveLineages;
import dr.evomodel.coalescent.StructuredCoalescentSchedule;
import dr.evomodel.treedatalikelihood.TreeTraversal;

import java.util.ArrayList;
import java.util.List;

import static dr.evomodel.coalescent.basta.ProcessOnCoalescentIntervalDelegate.*;

/**
 * @author Marc A Suchard
 * @author Guy Baele
 */
public class CoalescentIntervalTraversal extends TreeTraversal {


    public static class DegenerateIntervalException extends RuntimeException {
        public DegenerateIntervalException(String message) {
            super(message);
        }
    }


    public static class ActiveSetException extends RuntimeException {
        public final NodeRef node;
        public ActiveSetException(String message, NodeRef node) {
            super(message);
            this.node = node;
        }
    }

    private final BigFastTreeIntervals treeIntervals;
    private final BastaBranchIntervalOperationBuilder operationBuilder;

    private int currentLikelihoodInterval;

    private boolean slabMetadataRecording = false;
    private final BastaSlabMetadataBuilder slabMetadataBuilder = new BastaSlabMetadataBuilder();
    private static final int SLAB_DUMP_OPS_PER_BLOCK_DEFAULT = readSlabOpsPerBlockEnv();
    private static final int SLAB_DUMP_K_INDEX_OFFSET_PAT_DEFAULT = readKIndexOffsetPatEnv();
    private int slabOpsPerBlock = SLAB_DUMP_OPS_PER_BLOCK_DEFAULT;

    private int kIndexOffsetPat = SLAB_DUMP_K_INDEX_OFFSET_PAT_DEFAULT;

    private static int readSlabOpsPerBlockEnv() {
        try {
            String v = System.getenv("BASTA_SLAB_OPS_PER_BLOCK");
            if (v != null && !v.isEmpty()) {
                int parsed = Integer.parseInt(v.trim());
                if (parsed > 0) return parsed;
            }
        } catch (SecurityException | NumberFormatException ignored) {
            // Fall through to default.
        }
        return 32;
    }

    private static int readKIndexOffsetPatEnv() {
        try {
            String v = System.getenv("BASTA_K_INDEX_OFFSET_PAT");
            if (v != null && !v.isEmpty()) {
                int parsed = Integer.parseInt(v.trim());
                if (parsed > 0) return parsed;
            }
        } catch (SecurityException | NumberFormatException ignored) {
        }
        return 16;
    }

    protected CoalescentIntervalTraversal(final Tree tree,
                                          final BigFastTreeIntervals treeIntervals,
                                          final BranchRateModel branchRateModel,
                                          final int numberSubIntervals) {
        this(tree, treeIntervals, branchRateModel, numberSubIntervals, true);
    }

    protected CoalescentIntervalTraversal(final Tree tree,
                                          final BigFastTreeIntervals treeIntervals,
                                          final BranchRateModel branchRateModel,
                                          final int numberSubIntervals,
                                          final boolean checkForZeroLengthIntervals) {
        super(tree, branchRateModel, TraversalType.REVERSE_LEVEL_ORDER);

        this.treeIntervals = treeIntervals;
        this.operationBuilder = new BastaBranchIntervalOperationBuilder(
                tree, branchRateModel, numberSubIntervals, checkForZeroLengthIntervals,
                branchIntervalOperations, matrixOperations, intervalStarts);
    }

    @Override
    public final void dispatchTreeTraversalCollectBranchAndNodeOperations() {
        dispatchVersion++;
        matrixOperations.clear();
        branchIntervalOperations.clear();
        intervalStarts.clear();

        if (slabMetadataRecording) {
            slabMetadataBuilder.beginTraversal(
                    treeModel.getNodeCount(),
                    treeModel.getExternalNodeCount(),
                    treeModel.getExternalNodeCount(),
                    kIndexOffsetPat);
        }

        if (traversalType == TraversalType.REVERSE_LEVEL_ORDER) {
            traverseReverseCoalescentLevelOrder();
        } else {
            assert false : "Unknown traversal type";
        }
    }

    public void setSlabMetadataRecording(boolean enabled) {
        this.slabMetadataRecording = enabled;
    }

    public boolean isSlabMetadataRecording() {
        return slabMetadataRecording;
    }

    public int bufferForNode(int nodeNumber) {
        if (!slabMetadataRecording) {
            return nodeNumber;
        }
        return operationBuilder.bufferForNode(nodeNumber);
    }

    public int[] buildAndPackSlabMetadata(int slabOpsPerBlock) {
        if (!slabMetadataRecording) {
            throw new IllegalStateException(
                    "Slab metadata recording was not enabled.");
        }
        return slabMetadataBuilder.buildAndPackDirect(slabOpsPerBlock);
    }

    public int getLastPackedSlabMetadataLength() {
        return slabMetadataBuilder.getLastPackedLength();
    }

    public int getSlabOpsPerBlock() {
        return slabOpsPerBlock;
    }

    public int getKIndexOffsetPat() {
        return kIndexOffsetPat;
    }

    public void setSlabConstants(int slabOpsPerBlock, int kIndexOffsetPat) {
        if (slabOpsPerBlock > 0) {
            this.slabOpsPerBlock = slabOpsPerBlock;
        }
        if (kIndexOffsetPat > 0) {
            this.kIndexOffsetPat = kIndexOffsetPat;
        }
    }

    public List<BranchIntervalOperation> getBranchIntervalOperations() {
        return branchIntervalOperations;
    }

    public List<TransitionMatrixOperation> getMatrixOperations() {
        return matrixOperations;
    }

    public int getCoalescentIntervalCount() {
        return currentLikelihoodInterval + 1;
    }

    public List<Integer> getIntervalStarts() {
        return intervalStarts;
    }

    private void traverseReverseCoalescentLevelOrder() {
        // Zero-length coalescent intervals are rejected by the operation builder, as a DegenerateIntervalException
        final StructuredCoalescentSchedule schedule;
        try {
            schedule = StructuredCoalescentSchedule.fromTreeIntervals(treeModel, treeIntervals, true, false);
        } catch (StructuredCoalescentActiveLineages.InactiveLineageException e) {
            throw new ActiveSetException("Not in active set", treeModel.getNode(e.node));
        }
        operationBuilder.build(schedule, slabMetadataRecording ? slabMetadataBuilder : null);
        currentLikelihoodInterval = operationBuilder.getCoalescentIntervalCount() - 1;
    }

    private final List<BranchIntervalOperation> branchIntervalOperations = new ArrayList<>();
    private final List<TransitionMatrixOperation> matrixOperations = new ArrayList<>();
    private final List<Integer> intervalStarts = new ArrayList<>();

    private long dispatchVersion = 0L;

    public long getDispatchVersion() {
        return dispatchVersion;
    }
}
