/*
 * PatternListGenerator.java
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

package dr.app.beauti.generator;

import dr.app.beauti.components.ComponentFactory;
import dr.app.beauti.components.ancestralstates.AncestralStatesComponentOptions;
import dr.app.beauti.components.sequenceerror.SequenceErrorModelComponentOptions;
import dr.app.beauti.options.BeautiOptions;
import dr.app.beauti.options.PartitionData;
import dr.app.beauti.options.PartitionSubstitutionModel;
import dr.app.beauti.types.BinaryModelType;
import dr.app.beauti.util.XMLWriter;
import dr.evolution.alignment.Alignment;
import dr.evolution.alignment.SitePatterns;
import dr.evolution.datatype.DataType;
import dr.evoxml.*;
import dr.util.Attribute;
import dr.xml.XMLParser;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @author Walter Xie
 */
public class PatternListGenerator extends Generator {

    public PatternListGenerator(BeautiOptions options, ComponentFactory[] components) {
        super(options, components);
    }

    /**
     * Writes the pattern lists
     *
     * @param partition the partition data to write the pattern lists for
     * @param writer    the writer
     */
    public void writePatternList(PartitionData partition, XMLWriter writer) {
        writer.writeText("");

        AncestralStatesComponentOptions ancestralStatesOptions = (AncestralStatesComponentOptions) options
                .getComponentOptions(AncestralStatesComponentOptions.class);

        SequenceErrorModelComponentOptions sequenceErrorOptions = (SequenceErrorModelComponentOptions) options
                .getComponentOptions(SequenceErrorModelComponentOptions.class);

        PartitionSubstitutionModel model = partition.getPartitionSubstitutionModel();

        String codonHeteroPattern = model.getCodonHeteroPattern();
        int partitionCount = model.getCodonPartitionCount();

        boolean isAncestralStatesModel =  (!ancestralStatesOptions.usingAncestralStates(partition) &&
                !sequenceErrorOptions.usingSequenceErrorModel(partition));
        boolean isCovarionModel = model.getDataType().getType() == DataType.COVARION
                && model.getBinarySubstitutionModel() == BinaryModelType.BIN_COVARION;

        boolean unique = isAncestralStatesModel || isCovarionModel;
        boolean strip = isAncestralStatesModel || isCovarionModel;

        if (model.getDataType().getType() == DataType.NUCLEOTIDES && codonHeteroPattern != null && partitionCount > 1) {

            if (codonHeteroPattern.equals("112")) {
                writer.writeComment("The " + (unique ? "unique " : "") + "patterns for codon positions 1 & 2");
                writer.writeOpenTag(MergePatternsParser.MERGE_PATTERNS,
                        new Attribute[]{
                                new Attribute.Default<String>(XMLParser.ID, partition.getPrefix() + model.getPrefixCodon(1) + SitePatternsParser.PATTERNS),
                        }
                );
                writePatternList(partition, 0, 3, null, unique, strip, writer);
                writePatternList(partition, 1, 3, null, unique, strip, writer);

                writer.writeCloseTag(MergePatternsParser.MERGE_PATTERNS);

                writer.writeComment("The " + (unique ? "unique " : "") + "patterns for codon position 3");
                writePatternList(partition, 2, 3, model.getPrefixCodon(2), unique, strip, writer);

            } else {

                // pattern is 123
                for (int i = 1; i <= 3; i++) {
                    writer.writeComment("The " + (unique ? "unique " : "") + "patterns for codon position " + i);
                    writePatternList(partition, i - 1, 3, model.getPrefixCodon(i), unique, strip, writer);
                }

            }// END: pattern is 123

        } else {
            writePatternList(partition, 0, 1, "", unique, strip, writer);
        }
    }

    /**
     * Write a single pattern list
     *
     * @param partition the partition to write a pattern list for
     * @param offset    offset by
     * @param every     skip every
     * @param writer    the writer
     */
    private void writePatternList(final PartitionData partition, int offset, int every, final String codonPrefix, final boolean unique, final boolean strip, final XMLWriter writer) {

        Alignment alignment = partition.getAlignment();
        int from = partition.getFromSite();
        int to = partition.getToSite();
        int partEvery = partition.getEvery();
        if (partEvery > 1 && every > 1) throw new IllegalArgumentException();

        if (from < 1) from = 1;
        every = Math.max(partEvery, every);

        from += offset;

        // this object is created solely to calculate the number of patterns in the alignment
        SitePatterns patterns = new SitePatterns(alignment, null, from - 1, to - 1, every, strip,
                unique ? SitePatterns.CompressionType.UNIQUE_ONLY : SitePatterns.CompressionType.UNCOMPRESSED);

        writer.writeComment("The " + (unique ? "unique " : "") + "patterns from " + from + " to " + (to > 0 ? to : "end") + ((every > 1) ? " every " + every : ""),
                "npatterns=" + patterns.getPatternCount());

        List<Attribute> attributes = new ArrayList<Attribute>();

        // no 11 of 112 codon, which uses mergePatterns id instead
        if (codonPrefix != null) {
            attributes.add(new Attribute.Default<String>(XMLParser.ID, partition.getPrefix() + codonPrefix + SitePatternsParser.PATTERNS));
        }

        attributes.add(new Attribute.Default<String>("from", "" + from));
        if (to >= 0) attributes.add(new Attribute.Default<String>("to", "" + to));

        if (every > 1) {
            attributes.add(new Attribute.Default<String>("every", "" + every));
        }

        if(!unique) {
            attributes.add(new Attribute.Default<Boolean>(SitePatternsParser.UNIQUE, false));
        }

        if (strip) {
            attributes.add(new Attribute.Default<Boolean>(SitePatternsParser.STRIP, false)); // default true
        }

        // generate <patterns>
        writer.writeOpenTag(SitePatternsParser.PATTERNS, attributes);
        writer.writeIDref(AlignmentParser.ALIGNMENT, alignment.getId());
        writer.writeCloseTag(SitePatternsParser.PATTERNS);
    }

}
