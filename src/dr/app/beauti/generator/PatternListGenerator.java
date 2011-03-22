package dr.app.beauti.generator;

import dr.app.beauti.components.ComponentFactory;
import dr.app.beauti.options.BeautiOptions;
import dr.app.beauti.options.PartitionData;
import dr.app.beauti.options.PartitionPattern;
import dr.app.beauti.options.PartitionSubstitutionModel;
import dr.app.beauti.util.XMLWriter;
import dr.evolution.alignment.Alignment;
import dr.evolution.alignment.SitePatterns;
import dr.evolution.datatype.Microsatellite;
import dr.evolution.datatype.Nucleotides;
import dr.evoxml.*;
import dr.util.Attribute;
import dr.xml.XMLParser;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Alexei Drummond
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

        PartitionSubstitutionModel model = partition.getPartitionSubstitutionModel();

        String codonHeteroPattern = model.getCodonHeteroPattern();
        int partitionCount = model.getCodonPartitionCount();

        if (model.getDataType() == Nucleotides.INSTANCE && codonHeteroPattern != null && partitionCount > 1) {

            if (codonHeteroPattern.equals("112")) {
                writer.writeComment("The unique patterns for codon positions 1 & 2");
                writer.writeOpenTag(MergePatternsParser.MERGE_PATTERNS,
                        new Attribute[]{
                                new Attribute.Default<String>(XMLParser.ID, model.getPrefix(1) + partition.getPrefix() + SitePatternsParser.PATTERNS),
                        }
                );
                writePatternList(partition, 0, 3, null, writer);
                writePatternList(partition, 1, 3, null, writer);

                writer.writeCloseTag(MergePatternsParser.MERGE_PATTERNS);

                writer.writeComment("The unique patterns for codon positions 3");
//                writer.writeOpenTag(MergePatternsParser.MERGE_PATTERNS,
//                        new Attribute[]{
//                                new Attribute.Default<String>(XMLParser.ID, model.getPrefix(2) + partition.getPrefix() + SitePatternsParser.PATTERNS),
//                        }
//                );

                writePatternList(partition, 2, 3, model.getPrefix(2), writer);

//                writer.writeCloseTag(MergePatternsParser.MERGE_PATTERNS);

            } else {
                // pattern is 123
                // write pattern lists for all three codon positions
                for (int i = 1; i <= 3; i++) {
                    writer.writeComment("The unique patterns for codon positions " + i);
//                    writer.writeOpenTag(MergePatternsParser.MERGE_PATTERNS,
//                            new Attribute[]{
//                                    new Attribute.Default<String>(XMLParser.ID, model.getPrefix(i) + partition.getPrefix() + SitePatternsParser.PATTERNS),
//                            }
//                    );

                    writePatternList(partition, i - 1, 3, model.getPrefix(i), writer);

//                    writer.writeCloseTag(MergePatternsParser.MERGE_PATTERNS);
                }

            }
        } else {
            writePatternList(partition, 0, 1, null, writer);
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
    private void writePatternList(PartitionData partition, int offset, int every, String codonPrefix, XMLWriter writer) {

        Alignment alignment = partition.getAlignment();
        int from = partition.getFromSite();
        int to = partition.getToSite();
        int partEvery = partition.getEvery();
        if (partEvery > 1 && every > 1) throw new IllegalArgumentException();

        if (from < 1) from = 1;
        every = Math.max(partEvery, every);

        from += offset;


        // this object is created solely to calculate the number of patterns in the alignment
        SitePatterns patterns = new SitePatterns(alignment, from - 1, to - 1, every);

        writer.writeComment("The unique patterns from " + from + " to " + (to > 0 ? to : "end") + ((every > 1) ? " every " + every : ""),
                "npatterns=" + patterns.getPatternCount());

        List<Attribute> attributes = new ArrayList<Attribute>();

        // no codon, unique patterns site patterns
        if ((offset == 0 && every == 1) || (codonPrefix != null) )
            attributes.add(new Attribute.Default<String>(XMLParser.ID, codonPrefix + partition.getPrefix() + SitePatternsParser.PATTERNS));

        attributes.add(new Attribute.Default<String>("from", "" + from));
        if (to >= 0) attributes.add(new Attribute.Default<String>("to", "" + to));

        if (every > 1) {
            attributes.add(new Attribute.Default<String>("every", "" + every));
        }

        // generate <patterns>
        writer.writeOpenTag(SitePatternsParser.PATTERNS, attributes);
        writer.writeIDref(AlignmentParser.ALIGNMENT, alignment.getId());
        writer.writeCloseTag(SitePatternsParser.PATTERNS);
    }

    /**
     * Micro-sat
     * @param partition
     * @param microsatList
     * @param writer
     */
    public void writePatternList(PartitionPattern partition, List<Microsatellite> microsatList, XMLWriter writer) {

        PartitionSubstitutionModel model = partition.getPartitionSubstitutionModel();

        if (model.getDataType().getType() == Microsatellite.INSTANCE.getType()) {
            writer.writeComment("The patterns for microsatellite");
            writer.writeOpenTag(MicrosatellitePatternParser.MICROSATPATTERN,
                    new Attribute[]{
                            new Attribute.Default<String>(XMLParser.ID, partition.getName()),
                    });

            Microsatellite m = model.getMicrosatellite();

            if (!microsatList.contains(m)) {
                microsatList.add(m);
                writer.writeTag(MicrosatelliteParser.MICROSAT,
                        new Attribute[]{
                                new Attribute.Default<String>(XMLParser.ID, m.getName()),
                                new Attribute.Default<Integer>(MicrosatelliteParser.MAX, m.getMax()),
                                new Attribute.Default<Integer>(MicrosatelliteParser.MIN, m.getMin()),
                        }, true);

            } else {
                writer.writeTag(MicrosatelliteParser.MICROSAT,
                        new Attribute[]{
                                new Attribute.Default<String>(XMLParser.IDREF, m.getName()),
                        }, true);
            }

            writer.writeOpenTag(MicrosatellitePatternParser.MICROSAT_SEQ);
            String seq = "";
            for (int i = 0; i < partition.getPatterns().getPatternCount(); i++) {
                if (i > 0) seq +=",";
                seq += partition.getPatterns().getPattern(i);
            }
            writer.writeText(seq);
            writer.writeCloseTag(MicrosatellitePatternParser.MICROSATPATTERN);

            writer.writeCloseTag(MicrosatellitePatternParser.MICROSATPATTERN);

        } else {

        }
    }

}
