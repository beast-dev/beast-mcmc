package dr.app.beauti.generator;

import dr.app.beauti.components.ComponentFactory;
import dr.app.beauti.options.BeautiOptions;
import dr.app.beauti.options.PartitionTreeModel;
import dr.app.beauti.util.XMLWriter;
import dr.evolution.util.Taxa;
import dr.evomodel.tree.TreeModel;
import dr.evomodelxml.speciation.SpeciesTreeModelParser;
import dr.evomodelxml.tree.MonophylyStatisticParser;
import dr.evomodelxml.tree.TMRCAStatisticParser;
import dr.evoxml.TaxaParser;
import dr.evoxml.TaxonParser;
import dr.util.Attribute;
import dr.xml.XMLParser;

import java.util.List;
import java.util.Map;


/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @author Walter Xie
 */
public class TMRCAStatisticsGenerator extends Generator {


    public TMRCAStatisticsGenerator(BeautiOptions options, ComponentFactory[] components) {
        super(options, components);
    }

    /**
     * Generate additional taxon sets
     *
     * @param writer    the writer
     * @param taxonSets a list of taxa to write
     */
    public void writeTaxonSets(XMLWriter writer, List<Taxa> taxonSets) {
        writer.writeText("");

        for (Taxa taxa : taxonSets) {
            writer.writeOpenTag(
                    TaxaParser.TAXA,
                    new Attribute[]{
                            new Attribute.Default<String>(XMLParser.ID, taxa.getId())
                    }
            );

            for (int j = 0; j < taxa.getTaxonCount(); j++) {
                writer.writeIDref(TaxonParser.TAXON, taxa.getTaxon(j).getId());
            }
            writer.writeCloseTag(TaxaParser.TAXA);
        }
    }

    /**
     * Generate tmrca statistics
     *
     * @param writer       the writer
     */
    public void writeTMRCAStatistics(XMLWriter writer) {
        List<Taxa> taxonSets;
        Map<Taxa, Boolean> taxonSetsMono;

        if (options.useStarBEAST) {
            taxonSets = options.speciesSets;
            taxonSetsMono = options.speciesSetsMono;

            writer.writeComment("Species Sets");
            writer.writeText("");
            for (Taxa taxa : taxonSets) {
                writer.writeOpenTag(TMRCAStatisticParser.TMRCA_STATISTIC,
                        new Attribute[]{
                                new Attribute.Default<String>(XMLParser.ID, "tmrca(" + taxa.getId() + ")"),
//                        new Attribute.Default<Boolean>(TMRCAStatisticParser.STEM, options.taxonSetsIncludeStem.get(taxa)),
                        }
                ); // make tmrca(tree.name) eay to read in log for Tracer
                writer.writeOpenTag(TMRCAStatisticParser.MRCA);
                writer.writeIDref(TaxaParser.TAXA, taxa.getId());
                writer.writeCloseTag(TMRCAStatisticParser.MRCA);
                writer.writeIDref(SpeciesTreeModelParser.SPECIES_TREE, SP_TREE);
                writer.writeCloseTag(TMRCAStatisticParser.TMRCA_STATISTIC);

                if (taxonSetsMono.get(taxa)) {
//                    && treeModel.getPartitionTreePrior().getNodeHeightPrior() != TreePriorType.YULE
//                    && options.getKeysFromValue(options.taxonSetsTreeModel, treeModel).size() > 1) {
                    writer.writeOpenTag(
                            MonophylyStatisticParser.MONOPHYLY_STATISTIC,
                            new Attribute[]{
                                    new Attribute.Default<String>(XMLParser.ID, "monophyly(" + taxa.getId() + ")"),
                            });
                    writer.writeOpenTag(MonophylyStatisticParser.MRCA);
                    writer.writeIDref(TaxaParser.TAXA, taxa.getId());
                    writer.writeCloseTag(MonophylyStatisticParser.MRCA);
                    writer.writeIDref(SpeciesTreeModelParser.SPECIES_TREE, SP_TREE);
                    writer.writeCloseTag(MonophylyStatisticParser.MONOPHYLY_STATISTIC);
                }
            }

        } else {
            taxonSets = options.taxonSets;
            taxonSetsMono = options.taxonSetsMono;

            writer.writeComment("Taxon Sets");
            writer.writeText("");
            for (Taxa taxa : taxonSets) {
                PartitionTreeModel treeModel = options.taxonSetsTreeModel.get(taxa);
                writer.writeOpenTag(TMRCAStatisticParser.TMRCA_STATISTIC,
                        new Attribute[]{
                                new Attribute.Default<String>(XMLParser.ID, "tmrca(" + treeModel.getPrefix() + taxa.getId() + ")"),
                                new Attribute.Default<Boolean>(TMRCAStatisticParser.STEM, options.taxonSetsIncludeStem.get(taxa)),
                        }
                ); // make tmrca(tree.name) eay to read in log for Tracer
                writer.writeOpenTag(TMRCAStatisticParser.MRCA);
                writer.writeIDref(TaxaParser.TAXA, taxa.getId());
                writer.writeCloseTag(TMRCAStatisticParser.MRCA);
                writer.writeIDref(TreeModel.TREE_MODEL, treeModel.getPrefix() + TreeModel.TREE_MODEL);
                writer.writeCloseTag(TMRCAStatisticParser.TMRCA_STATISTIC);

                if (taxonSetsMono.get(taxa)) {
//                    && treeModel.getPartitionTreePrior().getNodeHeightPrior() != TreePriorType.YULE
//                    && options.getKeysFromValue(options.taxonSetsTreeModel, treeModel).size() > 1) {
                    writer.writeOpenTag(
                            MonophylyStatisticParser.MONOPHYLY_STATISTIC,
                            new Attribute[]{
                                    new Attribute.Default<String>(XMLParser.ID, "monophyly(" + taxa.getId() + ")"),
                            });
                    writer.writeOpenTag(MonophylyStatisticParser.MRCA);
                    writer.writeIDref(TaxaParser.TAXA, taxa.getId());
                    writer.writeCloseTag(MonophylyStatisticParser.MRCA);
                    writer.writeIDref(TreeModel.TREE_MODEL, treeModel.getPrefix() + TreeModel.TREE_MODEL);
                    writer.writeCloseTag(MonophylyStatisticParser.MONOPHYLY_STATISTIC);
                }
            }
        }
    }


}