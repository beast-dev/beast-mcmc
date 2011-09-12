package dr.app.beauti.generator;

import dr.app.beauti.components.ComponentFactory;
import dr.app.beauti.components.dnds.DnDsComponentOptions;
import dr.app.beauti.options.BeautiOptions;
import dr.app.beauti.options.PartitionData;
import dr.app.beauti.options.PartitionPattern;
import dr.app.beauti.options.PartitionSubstitutionModel;
import dr.app.beauti.types.BinaryModelType;
import dr.app.beauti.util.XMLWriter;
import dr.evolution.alignment.Alignment;
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

        boolean isCovarionModel = model.getDataType().getType() == DataType.COVARION
                && model.getBinarySubstitutionModel() == BinaryModelType.BIN_COVARION;

        if (model.getDataType() == Nucleotides.INSTANCE && codonHeteroPattern != null && partitionCount > 1) {

            if (codonHeteroPattern.equals("112")) {
                writer.writeComment("The unique patterns for codon positions 1 & 2");
                writer.writeOpenTag(MergePatternsParser.MERGE_PATTERNS,
                        new Attribute[]{
                                new Attribute.Default<String>(XMLParser.ID, model.getPrefix(1) + partition.getPrefix() + SitePatternsParser.PATTERNS),
                        }
                );
                writePatternList(partition, 0, 3, null, isCovarionModel, writer);
                writePatternList(partition, 1, 3, null, isCovarionModel, writer);

                writer.writeCloseTag(MergePatternsParser.MERGE_PATTERNS);

                writer.writeComment("The unique patterns for codon positions 3");
//                writer.writeOpenTag(MergePatternsParser.MERGE_PATTERNS,
//                        new Attribute[]{
//                                new Attribute.Default<String>(XMLParser.ID, model.getPrefix(2) + partition.getPrefix() + SitePatternsParser.PATTERNS),
//                        }
//                );

                writePatternList(partition, 2, 3, model.getPrefix(2), isCovarionModel, writer);

//                writer.writeCloseTag(MergePatternsParser.MERGE_PATTERNS);

            } else {

        		DnDsComponentOptions component = (DnDsComponentOptions) options
				.getComponentOptions(DnDsComponentOptions.class);

				boolean doRobustCounting = component.doRobustCounting();

				if (doRobustCounting) {

//					System.out.println("HERE");
					for (int i = 1; i <= 3; i++) {
						writePatternList(partition, i - 1, 3, model
								.getPrefix(i), false, writer);
					}

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

                    writePatternList(partition, i - 1, 3, model.getPrefix(i), isCovarionModel, writer);

						// writer.writeCloseTag(MergePatternsParser.MERGE_PATTERNS);
					}
				}// END: doRobustCounting

			}// END: pattern is 123

        } else {
            writePatternList(partition, 0, 1, "", isCovarionModel, writer);
        }
    }

    private void writePatternList(PartitionData partition, int offset, int every, String codonPrefix, boolean isCovarionModel, XMLWriter writer) {
    	writePatternList(partition, offset, every, codonPrefix, true, isCovarionModel, writer);
    }
    /**
     * Write a single pattern list
     *
     * @param partition the partition to write a pattern list for
     * @param offset    offset by
     * @param every     skip every
     * @param writer    the writer
     */
    private void writePatternList(PartitionData partition, int offset, int every, String codonPrefix, boolean unique, boolean isCovarionModel, XMLWriter writer) {

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

        // no 11 of 112 codon, which uses mergePatterns
        if (codonPrefix != null)
            attributes.add(new Attribute.Default<String>(XMLParser.ID, codonPrefix + partition.getPrefix() + SitePatternsParser.PATTERNS));

        attributes.add(new Attribute.Default<String>("from", "" + from));
        if (to >= 0) attributes.add(new Attribute.Default<String>("to", "" + to));

        if (every > 1) {
            attributes.add(new Attribute.Default<String>("every", "" + every));
        }

        if(!unique) {
        	attributes.add(new Attribute.Default<Boolean>(SitePatternsParser.UNIQUE, false));
        }

        if (isCovarionModel) {
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
            writer.writeComment("The patterns for microsatellite");
            writer.writeOpenTag(MicrosatellitePatternParser.MICROSATPATTERN,
                    new Attribute[]{
                            new Attribute.Default<String>(XMLParser.ID, partition.getName()),
                    });

            if (options.partitionsHaveIdenticalTaxa()) {
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
                                new Attribute.Default<Integer>(MicrosatelliteParser.MICROSAT, m.getMin()),
                        }, true);

            } else {
                writer.writeTag(MicrosatelliteParser.MICROSAT,
                        new Attribute[]{
                                new Attribute.Default<String>(XMLParser.IDREF, m.getName()),
                        }, true);
            }

            writer.writeOpenTag(MicrosatellitePatternParser.MICROSAT_SEQ);
            String seq = "";
            for (int i = 0; i < partition.getPatterns().getTaxonCount(); i++) {
                if (i > 0) seq +=",";
                seq += partition.getPatterns().getPatternState(i, 0);
            }
            writer.writeText(seq);
            writer.writeCloseTag(MicrosatellitePatternParser.MICROSAT_SEQ);

            writer.writeCloseTag(MicrosatellitePatternParser.MICROSATPATTERN);

        } else {

        }
    }

}
