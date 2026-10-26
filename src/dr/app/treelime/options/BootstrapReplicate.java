/*
 * BootstrapReplicate.java
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

package dr.app.treelime.options;

import dr.evolution.alignment.Alignment;

import java.util.HashMap;
import java.util.Map;

/**
 * One bootstrap replicate's data: per partition, an alignment of its resampled columns (read whole by the
 * partition's patterns) and, for a partition padded with constant sites, the redrawn constant-site counts. The
 * generator writes these alignments next to the originals, which the starting tree still uses, so every replicate
 * keeps the original fit's topology.
 */
public class BootstrapReplicate {

    private final Map<PartitionData, Alignment> alignments = new HashMap<PartitionData, Alignment>();
    private final Map<PartitionData, double[]> constantSiteCounts = new HashMap<PartitionData, double[]>();

    public void put(PartitionData partition, Alignment alignment, double[] constantCounts) {
        alignments.put(partition, alignment);
        if (constantCounts != null) {
            constantSiteCounts.put(partition, constantCounts.clone());
        }
    }

    public Alignment getAlignment(PartitionData partition) {
        return alignments.get(partition);
    }

    /**
     * @return the redrawn constant-site counts, or null when the partition is not padded
     */
    public double[] getConstantSiteCounts(PartitionData partition) {
        return constantSiteCounts.get(partition);
    }

    /**
     * @return the id of the partition's resampled alignment in the XML
     */
    public static String alignmentId(PartitionData partition) {
        return partition.getPrefix() + "bootstrap.alignment";
    }
}
