package dr.app.beauti.generator;

import dr.app.beauti.components.ComponentFactory;
import dr.app.beauti.enumTypes.BinaryModelType;
import dr.app.beauti.options.BeautiOptions;
import dr.app.beauti.options.PartitionData;
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

import java.util.ArrayList;
import java.util.List;


/**
 * @author Alexei Drummond
 * @author Walter Xie
 */
public class AlignmentGenerator extends Generator {


    public AlignmentGenerator(BeautiOptions options, ComponentFactory[] components) {
        super(options, components);
    }

    /**
     * Write all Alignments
     * @param writer  XMLWriter
     */
    public void writeAlignments(XMLWriter writer) {
        List<Alignment> alignments = new ArrayList<Alignment>();

        for (PartitionData partition : options.getNonTraitsDataList()) {
            Alignment alignment = partition.getAlignment();
            if (!alignments.contains(alignment)) {
                alignments.add(alignment);
            }
        }

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
        writer.writeComment("The sequence alignment (each sequence refers to a taxon above).");
        writer.writeComment("ntax=" + alignment.getTaxonCount() + " nchar=" + alignment.getSiteCount());
        if (options.samplePriorOnly) {
            writer.writeComment("Null sequences generated in order to sample from the prior only.");
        }

        writer.writeOpenTag(
                AlignmentParser.ALIGNMENT,
                new Attribute[]{
                        new Attribute.Default<String>(XMLParser.ID, alignment.getId()),
                        new Attribute.Default<String>("dataType", getAlignmentDataTypeDescription(alignment))
                }
        );

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
                writer.writeText("N");
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
                }
                break;
        }

        return description;
    }

}