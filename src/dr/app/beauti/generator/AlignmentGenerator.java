/*
 * AlignmentGenerator.java
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

package dr.app.beauti.generator;

import dr.app.beauti.components.ComponentFactory;
import dr.app.beauti.components.dollo.DolloComponentOptions;
import dr.app.beauti.options.BeautiOptions;
import dr.app.beauti.types.BinaryModelType;
import dr.app.beauti.util.XMLWriter;
import dr.evolution.alignment.Alignment;
import dr.evolution.datatype.DataType;
import dr.evolution.datatype.TwoStateCovarion;
import dr.evolution.util.Taxon;
import dr.evoxml.AlignmentParser;
import dr.evoxml.SequenceParser;
import dr.evoxml.TaxonParser;
import dr.util.Attribute;
import dr.xml.XMLParser;

import java.util.List;


/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @author Walter Xie
 */
public class AlignmentGenerator extends Generator {

    // used as a null sequence if sampling from the prior only. It is
    // three characters long in case it is codon partitioned.
    private static final String NULL_SEQUENCE = "???";

    public AlignmentGenerator(BeautiOptions options, ComponentFactory[] components) {
        super(options, components);
    }

    /**
     * Write all Alignments
     * @param writer  XMLWriter
     */
    public void writeAlignments(List<Alignment> alignments, XMLWriter writer) {
        int index = 1;
        for (Alignment alignment : alignments) {
            if (alignments.size() > 1) {
                //if (!options.allowDifferentTaxa) {
                alignment.setId(AlignmentParser.ALIGNMENT + index);
                //} else { // e.g. alignment_gene1
                // alignment.setId("alignment_" + mulitTaxaTagName + index);
                //}
            } else {
                alignment.setId(AlignmentParser.ALIGNMENT);
            }
            writeAlignment(alignment, writer);
            index += 1;
            writer.writeText("");
        }

    }

    /**
     * Generate an alignment block from these beast options
     *
     * @param alignment the alignment to write
     * @param writer    the writer
     */
    private void writeAlignment(Alignment alignment, XMLWriter writer) {

        writer.writeText("");
        writer.writeComment("The sequence alignment (each sequence refers to a taxon above).",
                "ntax=" + alignment.getTaxonCount() + " nchar=" + alignment.getSiteCount());
        if (options.samplePriorOnly) {
            writer.writeComment("Null sequences generated in order to sample from the prior only.");
        }

        if (getAlignmentDataTypeDescription(alignment) != null) {
            writer.writeOpenTag(
                    AlignmentParser.ALIGNMENT,
                    new Attribute[]{
                            new Attribute.Default<String>(XMLParser.ID, alignment.getId()),
                            new Attribute.Default<String>(DataType.DATA_TYPE, getAlignmentDataTypeDescription(alignment))
                    }
            );
        } else {
            writer.writeOpenTag(
                    AlignmentParser.ALIGNMENT, new Attribute.Default<String>(XMLParser.ID, alignment.getId()));

            writer.writeIDref(DataType.DATA_TYPE, getAlignmentDataTypeIdref(alignment));
        }

        for (int i = 0; i < alignment.getTaxonCount(); i++) {
            Taxon taxon = alignment.getTaxon(i);

            writer.writeOpenTag(SequenceParser.SEQUENCE);
            writer.writeIDref(TaxonParser.TAXON, taxon.getId());
            if (!options.samplePriorOnly) {
//                writer.checkText(alignment.getAlignedSequenceString(i));
                writer.writeText(alignment.getAlignedSequenceString(i));

//                System.out.println(taxon.getId() + ": \n" + alignment.getAlignedSequenceString(i));
//                System.out.println("len = " + alignment.getAlignedSequenceString(i).length() + "\n");
            } else {
                // generate a codon in case there is codon partitioning
                writer.writeText(NULL_SEQUENCE);
            }
            writer.writeCloseTag(SequenceParser.SEQUENCE);
        }
        writer.writeCloseTag(AlignmentParser.ALIGNMENT);
    }


    /**
     * Determine and return the datatype description for these beast options
     * note that the datatype in XML may differ from the actual datatype
     *
     * @param alignment the alignment to get data type description of
     * @return description
     */

    private String getAlignmentDataTypeDescription(Alignment alignment) {
        String description = alignment.getDataType().getDescription();

        switch (alignment.getDataType().getType()) {
            case DataType.TWO_STATES: // when dataType="binary"
//            case DataType.COVARION:
//                description = alignment.getDataType().getDescription();

                // if choose Covarion model then should change into dataType="twoStateCovarion"
                if (options.getPartitionData(alignment).getPartitionSubstitutionModel().getBinarySubstitutionModel()
                        == BinaryModelType.BIN_COVARION) {
                    description = TwoStateCovarion.INSTANCE.getDescription(); // dataType="twoStateCovarion"
                } else if (options.getPartitionData(alignment).getPartitionSubstitutionModel().getBinarySubstitutionModel()
                        == BinaryModelType.BIN_DOLLO) {
                    description = null;
                }
                break;
        }

        return description;
    }

    private String getAlignmentDataTypeIdref(Alignment alignment) {
        String description = alignment.getDataType().getDescription();
        switch (alignment.getDataType().getType()) {
            case DataType.TWO_STATES:
                if (options.getPartitionData(alignment).getPartitionSubstitutionModel().getBinarySubstitutionModel()
                        == BinaryModelType.BIN_DOLLO) {
                    description = DolloComponentOptions.DATA_NAME;
                }
        }
        return description;
    }
}