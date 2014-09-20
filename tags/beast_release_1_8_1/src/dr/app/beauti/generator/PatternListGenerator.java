package dr.app.beauti.generator;

import dr.app.beauti.components.ComponentFactory;
import dr.app.beauti.components.ancestralstates.AncestralStatesComponentOptions;
import dr.app.beauti.components.sequenceerror.SequenceErrorModelComponentOptions;
import dr.app.beauti.options.BeautiOptions;
import dr.app.beauti.options.PartitionData;
import dr.app.beauti.options.PartitionPattern;
import dr.app.beauti.options.PartitionSubstitutionModel;
import dr.app.beauti.types.BinaryModelType;
import dr.app.beauti.util.XMLWriter;
import dr.evolution.alignment.Alignment;
import dr.evolution.alignment.Patterns;
import dr.evolution.alignment.SitePatterns;
import dr.evolution.datatype.DataType;
import dr.evolution.datatype.Microsatellite;
import dr.evolution.datatype.Nucleotides;
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
        SitePatterns patterns = new SitePatterns(alignment, null, from - 1, to - 1, every, strip, unique);

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

    /**
     * Micro-sat
     * @param partition
     * @param microsatList
     * @param writer
     */
    public void writePatternList(PartitionPattern partition, List<Microsatellite> microsatList, XMLWriter writer) throws GeneratorException {

        PartitionSubstitutionModel model = partition.getPartitionSubstitutionModel();

        if (model.getDataType().getType() == DataType.MICRO_SAT) {
            Patterns patterns = partition.getPatterns();

            writer.writeComment("The patterns for microsatellite");
            writer.writeOpenTag(MicrosatellitePatternParser.MICROSATPATTERN,
                    new Attribute[]{
                            new Attribute.Default<String>(XMLParser.ID, partition.getName()),
                    });

            if (options.hasIdenticalTaxa() && !patterns.hasMask()) {
                writer.writeIDref(TaxaParser.TAXA, TaxaParser.TAXA);
            } else {
                writer.writeIDref(TaxaParser.TAXA, partition.getName() + "." + TaxaParser.TAXA);
            }

            Microsatellite m = model.getMicrosatellite();
            if (m == null) throw new GeneratorException("Microsatellite is null in partition:\n" + partition.getName());

            if (!microsatList.contains(m)) {
                microsatList.add(m);
                writer.writeTag(MicrosatelliteParser.MICROSAT,
                        new Attribute[]{
                                new Attribute.Default<String>(XMLParser.ID, m.getName()),
                                new Attribute.Default<Integer>(MicrosatelliteParser.MAX, m.getMax()),
                                new Attribute.Default<Integer>(MicrosatelliteParser.MIN, m.getMin()),
                                new Attribute.Default<Integer>(MicrosatelliteParser.UNIT_LENGTH, m.getUnitLength()),
                        }, true);

            } else {
                writer.writeTag(MicrosatelliteParser.MICROSAT,
                        new Attribute[]{
                                new Attribute.Default<String>(XMLParser.IDREF, m.getName()),
                        }, true);
            }

            writer.writeOpenTag(MicrosatellitePatternParser.MICROSAT_SEQ);

            String seq = "";
            int c = 0;
            for (int i = 0; i < patterns.getTaxonCount(); i++) {
                if (!patterns.isMasked(i)) {
                    if (c > 0) seq += ",";
                    int state = patterns.getPatternState(i, 0);
                    if (state == Microsatellite.UNKNOWN_STATE_LENGTH) {
                        seq += Microsatellite.UNKNOWN_CHARACTER;
                    } else {
                        seq += Integer.toString(state);
                    }
                    c++;
                }
            }
            writer.writeText(seq);

            writer.writeCloseTag(MicrosatellitePatternParser.MICROSAT_SEQ);

            writer.writeCloseTag(MicrosatellitePatternParser.MICROSATPATTERN);

        } else {

        }
    }

}
