package dr.app.beauti.generator;

import dr.app.beauti.components.ComponentFactory;
import dr.app.beauti.options.*;
import dr.app.beauti.util.XMLWriter;
import dr.evomodel.tree.TMRCAStatistic;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.tree.MonophylyStatistic;
import dr.evoxml.*;
import dr.util.Attribute;
import dr.xml.XMLParser;
import dr.evolution.util.Taxa;

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
            writer.writeOpenTag(TMRCAStatistic.TMRCA_STATISTIC,
                new Attribute[]{new Attribute.Default<String>(XMLParser.ID, "tmrca("
                        + taxa.getTreeModel().getPrefix() + taxa.getId() + ")"),}
            ); // make tmrca(tree.name) eay to read in log for Tracer
            writer.writeOpenTag(TMRCAStatistic.MRCA);
            writer.writeIDref(TaxaParser.TAXA, taxa.getId());
            writer.writeCloseTag(TMRCAStatistic.MRCA);
            writer.writeIDref(TreeModel.TREE_MODEL, taxa.getTreeModel().getPrefix() + TreeModel.TREE_MODEL);
            writer.writeCloseTag(TMRCAStatistic.TMRCA_STATISTIC);

            if (options.taxonSetsMono.get(taxa)) {
                writer.writeOpenTag(
                        MonophylyStatistic.MONOPHYLY_STATISTIC,
                        new Attribute[]{
                                new Attribute.Default<String>(XMLParser.ID, "monophyly(" + taxa.getId() + ")"),
                        });
                writer.writeOpenTag(MonophylyStatistic.MRCA);
                writer.writeIDref(TaxaParser.TAXA, taxa.getId());
                writer.writeCloseTag(MonophylyStatistic.MRCA);
                writer.writeIDref(TreeModel.TREE_MODEL, taxa.getTreeModel().getPrefix() + TreeModel.TREE_MODEL);
                writer.writeCloseTag(MonophylyStatistic.MONOPHYLY_STATISTIC);
            }
        }
    }


}