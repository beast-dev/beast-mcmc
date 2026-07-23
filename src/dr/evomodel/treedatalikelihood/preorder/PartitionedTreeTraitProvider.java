/*
 * PartitionedTreeTraitProvider.java
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

package dr.evomodel.treedatalikelihood.preorder;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evolution.tree.TreeTraitProvider;

public class PartitionedTreeTraitProvider implements TreeTraitProvider {

    private final TreeTrait[] originalTraits;
    private final TreeTrait[] partitionedTraits;
    private final int nPartitions;
    private final int[] partitionDimensions;
    private final int totalDim;

    private final int[] nodeOffsets;


    public PartitionedTreeTraitProvider(TreeTrait[] originalTraits, int[] partitionDimensions) {
        this.originalTraits = originalTraits;
        this.partitionDimensions = partitionDimensions;
        this.nPartitions = partitionDimensions.length;

        this.nodeOffsets = new int[nPartitions];
        int totalDim = 0;
        for (int i = 0; i < nPartitions; i++) {
            nodeOffsets[i] = totalDim;
            totalDim += partitionDimensions[i];
        }

        this.totalDim = totalDim;

        int nTraits = originalTraits.length;

        this.partitionedTraits = new TreeTrait[nTraits * nPartitions];


        for (int i = 0; i < nTraits; i++) {
            final int ind = i;

            for (int j = 0; j < nPartitions; j++) {
                final int partition = j;
                TreeTrait.DA trait = new TreeTrait.DA() {
                    @Override
                    public String getTraitName() {
                        return originalTraits[ind].getTraitName() + "." + partition; //TODO: could make this more descriptive
                    }

                    @Override
                    public Intent getIntent() {
                        return originalTraits[ind].getIntent();
                    }

                    @Override
                    public double[] getTrait(Tree tree, NodeRef node) {
                        double[] originalValues = (double[]) originalTraits[ind].getTrait(tree, node);
                        return partitionValues(originalValues, partition);
                    }
                };

                partitionedTraits[i * nPartitions + j] = trait;
            }
        }

    }

    private double[] partitionValues(double[] originalValues, int partition) {
        int nNodes = originalValues.length / totalDim;
        int partitionDim = partitionDimensions[partition];
        if (originalValues.length % totalDim != 0) {
            throw new RuntimeException("The trait dimension must be a factor of the tree trait length.");
        }

        double[] partitionedValues = new double[nNodes * partitionDimensions[partition]];
        int originalOffset = nodeOffsets[partition];
        int partitionedOffset = 0;
        for (int i = 0; i < nNodes; i++) {
            System.arraycopy(originalValues, originalOffset, partitionedValues, partitionedOffset, partitionDim);

            originalOffset += totalDim;
            partitionedOffset += partitionDim;
        }

        return partitionedValues;
    }


    @Override
    public TreeTrait[] getTreeTraits() {
        return partitionedTraits;
    }

    @Override
    public TreeTrait getTreeTrait(String key) {
        throw new RuntimeException("Error: not implemented.");
    }


}
