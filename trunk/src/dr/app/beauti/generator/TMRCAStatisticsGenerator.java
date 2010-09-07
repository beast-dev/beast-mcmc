package dr.app.beauti.generator;

import dr.app.beauti.components.ComponentFactory;
import dr.app.beauti.options.BeautiOptions;
import dr.app.beauti.util.XMLWriter;
import dr.evolution.util.Taxa;
import dr.evomodel.tree.TreeModel;
import dr.evomodelxml.tree.MonophylyStatisticParser;
import dr.evomodelxml.tree.TMRCAStatisticParser;
import dr.evoxml.TaxaParser;
import dr.evoxml.TaxonParser;
import dr.util.Attribute;
import dr.xml.XMLParser;

import java.util.List;


/**
 * @author Alexei Drummond
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
     * @param writer the writer
     */
    public void writeTMRCAStatistics(XMLWriter writer) {

        writer.writeText("");
        for (Taxa taxa : options.taxonSets) {
            writer.writeOpenTag(TMRCAStatisticParser.TMRCA_STATISTIC,
                new Attribute[]{
                    new Attribute.Default<String>(XMLParser.ID, "tmrca(" + taxa.getTreeModel().getPrefix() + taxa.getId() + ")"),
                    new Attribute.Default<Boolean>(TMRCAStatisticParser.PARENT, options.taxonSetsForParent.get(taxa)),
                }
            ); // make tmrca(tree.name) eay to read in log for Tracer
            writer.writeOpenTag(TMRCAStatisticParser.MRCA);
            writer.writeIDref(TaxaParser.TAXA, taxa.getId());
            writer.writeCloseTag(TMRCAStatisticParser.MRCA);
            writer.writeIDref(TreeModel.TREE_MODEL, taxa.getTreeModel().getPrefix() + TreeModel.TREE_MODEL);
            writer.writeCloseTag(TMRCAStatisticParser.TMRCA_STATISTIC);

            if (options.taxonSetsMono.get(taxa)) {
                writer.writeOpenTag(
                        MonophylyStatisticParser.MONOPHYLY_STATISTIC,
                        new Attribute[]{
                                new Attribute.Default<String>(XMLParser.ID, "monophyly(" + taxa.getId() + ")"),
                        });
                writer.writeOpenTag(MonophylyStatisticParser.MRCA);
                writer.writeIDref(TaxaParser.TAXA, taxa.getId());
                writer.writeCloseTag(MonophylyStatisticParser.MRCA);
                writer.writeIDref(TreeModel.TREE_MODEL, taxa.getTreeModel().getPrefix() + TreeModel.TREE_MODEL);
                writer.writeCloseTag(MonophylyStatisticParser.MONOPHYLY_STATISTIC);
            }
        }
    }


}