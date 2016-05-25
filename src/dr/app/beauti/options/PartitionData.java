/*
 * PartitionData.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
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
 */

package dr.app.beauti.options;

import dr.evolution.alignment.Alignment;
import dr.evolution.alignment.Patterns;
import dr.evolution.datatype.DataType;
import dr.evolution.util.TaxonList;

import java.util.List;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @author Walter Xie
 */
public class PartitionData extends AbstractPartitionData {

    private static final long serialVersionUID = 1642891822797102561L;

    private final Alignment alignment;

    private int fromSite;
    private int toSite;
    private int every = 1;


    public PartitionData(BeautiOptions options, String name, String fileName, Alignment alignment) {
        this(options, name, fileName, alignment, -1, -1, 1);
    }

    public PartitionData(BeautiOptions options, String name, String fileName, Alignment alignment, int fromSite, int toSite, int every) {
        super(options, name, fileName);
        this.alignment = alignment;

        this.fromSite = fromSite;
        this.toSite = toSite;
        this.every = every;

        this.traits = null;

        Patterns patterns = null;
        if (alignment != null) {
            patterns = new Patterns(alignment);
        }

        // This is too slow to be done at data loading.
        // calculateMeanDistance(patterns);
        calculateMeanDistance(null);
    }

    public PartitionData(BeautiOptions options, String name, String fileName, List<TraitData> traits) {
        super(options, name, fileName);
        this.alignment = null;

        this.fromSite = -1;
        this.toSite = -1;
        this.every = 1;

        this.traits = traits;

        calculateMeanDistance(null);
    }

    public Alignment getAlignment() {
        return alignment;
    }

    public int getFromSite() {
        return fromSite;
    }

    public int getToSite() {
        return toSite;
    }

    public int getEvery() {
        return every;
    }

    public TaxonList getTaxonList() {
        if (traits != null) {
//            // if this is a trait then just give the complete taxon list (taxa without specified
//            // traits are treated as missing data.
//            return options.taxonList;
            return null;
        }
        return getAlignment();
    }

    public int getSiteCount() {
        if (alignment != null) {
            int from = getFromSite();
            if (from < 1) {
                from = 1;
            }
            int to = getToSite();
            if (to < 1) {
                to = alignment.getSiteCount();
            }
            return (to - from + 1) / every;
        } else {
            return traits.size();
        }
    }

    public DataType getDataType() {
        if (alignment != null) {
            return alignment.getDataType();
        } else if (traits != null) {
            return traits.get(0).getDataType();
        } else {
            throw new RuntimeException("Trait and alignment are null");
        }
    }

    public String getDataDescription() {
        if (alignment != null) {
            return alignment.getDataType().getDescription();
        } else if (traits != null) {
            return traits.get(0).getTraitType().toString();
        } else {
            throw new RuntimeException("Trait and alignment are null");
        }
    }

    public String getPrefix() {
        String prefix = "";

        if (getTraits() != null) {
            // if it is a trait partition then always give a prefix
            prefix += getName() + ".";
        } else {
            // this method provides prefix as long as multi-data-partitions case,
            // because options.dataPartitions may contain traits, use options.getPartitionData()
            if (options.getPartitionData().size() > 1) { // getPartitionData() already excludes traits and microsatellite
                prefix += getName() + ".";
            }
        }

        return prefix;
    }

}
