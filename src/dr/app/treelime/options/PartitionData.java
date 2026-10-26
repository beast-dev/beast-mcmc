/*
 * PartitionData.java
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
import dr.evolution.alignment.PatternList;
import dr.evolution.alignment.Patterns;
import dr.evolution.alignment.SitePatterns;
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
    private PatternList patterns;
    private final int fromSite;
    private final int toSite;
    private final int every;
    // counts of the constant A, C, G and T columns left out of a SNP alignment (IQ-TREE -fconst); null for none
    private double[] constantSiteCounts = null;


    public PartitionData(TreeLimeOptions options, String name, String fileName, PatternList patterns) {
        super(options, name, fileName);
        this.alignment = null;
        this.patterns = patterns;

        this.fromSite = -1;
        this.toSite = -1;
        this.every = 1;

        this.traits = null;

        // This is too slow to be done at data loading.
        // calculateMeanDistance(patterns);
        calculateMeanDistance(null);
    }

    public PartitionData(TreeLimeOptions options, String name, String fileName, Alignment alignment, int fromSite, int toSite, int every) {
        super(options, name, fileName);
        this.alignment = alignment;

        this.fromSite = fromSite;
        this.toSite = toSite;
        this.every = every;

        this.traits = null;

        // This is too slow to be done at data loading.
        // calculateMeanDistance(patterns);
        calculateMeanDistance(null);
    }

    public PartitionData(TreeLimeOptions options, String name, String fileName, List<TraitData> traits) {
        super(options, name, fileName);
        this.alignment = null;

        this.fromSite = -1;
        this.toSite = -1;
        this.every = 1;

        this.traits = traits;

        calculateMeanDistance(null);
    }

    public PartitionData(TreeLimeOptions options, String name, String fileName, Alignment alignment) {
        this(options, name, fileName, alignment, -1, -1, 1);
    }


    public Alignment getAlignment() {
        return alignment;
    }

    /**
     * @return the counts of constant columns, one per state (A, C, G, T), added to the patterns as weights, or null
     */
    public double[] getConstantSiteCounts() {
        return constantSiteCounts;
    }

    /**
     * @param counts one count per state, or null for none
     */
    public void setConstantSiteCounts(double[] counts) {
        this.constantSiteCounts = counts == null ? null : counts.clone();
    }

    /**
     * Reads constant-site counts as typed or as written by snp-sites -C: four whole numbers for A, C, G and T,
     * separated by commas or white space.
     *
     * @return the counts, or null for blank text (no counts)
     * @throws IllegalArgumentException for anything else, with a message for the user
     */
    public static double[] parseConstantSiteCounts(String text) {
        String trimmed = text == null ? "" : text.trim();
        if (trimmed.isEmpty()) {
            return null;
        }
        String[] fields = trimmed.split("[,\\s]+");
        if (fields.length != 4) {
            throw new IllegalArgumentException("Expected four counts, for A, C, G and T, but found " + fields.length +
                    ": " + trimmed);
        }
        double[] counts = new double[4];
        for (int i = 0; i < 4; i++) {
            try {
                counts[i] = Double.parseDouble(fields[i]);
            } catch (NumberFormatException e) {
                throw new IllegalArgumentException("Not a count: " + fields[i]);
            }
            if (!(counts[i] >= 0) || Double.isInfinite(counts[i]) || counts[i] != Math.rint(counts[i])) {
                throw new IllegalArgumentException("Counts must be whole numbers, not negative: " + fields[i]);
            }
        }
        return counts;
    }

    /**
     * @return the number of this partition's columns where every taxon has the same, unambiguous state
     */
    public int getConstantColumnCount() {
        if (alignment == null) {
            return 0;
        }
        int from = fromSite < 1 ? 0 : fromSite - 1;
        int to = toSite < 1 ? alignment.getSiteCount() : toSite;
        int count = 0;
        for (int site = from; site < to; site += every) {
            int state = -1;
            boolean constant = true;
            for (int taxon = 0; taxon < alignment.getSequenceCount() && constant; taxon++) {
                int s = alignment.getState(taxon, site);
                if (alignment.getDataType().isAmbiguousState(s) || (state >= 0 && s != state)) {
                    constant = false;
                }
                state = s;
            }
            if (constant) {
                count++;
            }
        }
        return count;
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

    public void compressPatterns(SitePatterns.CompressionType compressionType) {
        assert alignment != null;

        patterns = new SitePatterns(alignment, alignment, compressionType);
    }


    public int getPatternCount() {
        if (alignment == null) {
            return traits.size();
        }
        if (patterns == null) {
            return -1;
//            patterns = new SitePatterns(alignment);
        }
        return patterns.getPatternCount();
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
            if (options.getPartitionData().size() > 1) { // getPartitionData() already excludes traits
                prefix += getName() + ".";
            }
        }

        return prefix;
    }

}
