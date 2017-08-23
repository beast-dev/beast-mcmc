/*
 * TMRCAStatisticsGenerator.java
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
                String id = "tmrca(" + treeModel.getPrefix() + taxa.getId() + ")";
                writeTMRCAStatistic(writer, id, taxa, treeModel, false, options.taxonSetsIncludeStem.get(taxa));

                if (treeModel.hasTipCalibrations()) {
                    id = "age(" + treeModel.getPrefix() + taxa.getId() + ")";
                    writeTMRCAStatistic(writer, id, taxa, treeModel, true, options.taxonSetsIncludeStem.get(taxa));
                }

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

    private void writeTMRCAStatistic(XMLWriter writer, String id, Taxa taxa, PartitionTreeModel treeModel, boolean isAbsolute, boolean includeStem) {
        writer.writeOpenTag(TMRCAStatisticParser.TMRCA_STATISTIC,
                new Attribute[]{
                        new Attribute.Default<String>(XMLParser.ID, id),
                        new Attribute.Default<Boolean>(TMRCAStatisticParser.ABSOLUTE, isAbsolute),
                        new Attribute.Default<Boolean>(TMRCAStatisticParser.STEM, includeStem),
                }
        ); // make tmrca(tree.name) eay to read in log for Tracer
        writer.writeOpenTag(TMRCAStatisticParser.MRCA);
        writer.writeIDref(TaxaParser.TAXA, taxa.getId());
        writer.writeCloseTag(TMRCAStatisticParser.MRCA);
        writer.writeIDref(TreeModel.TREE_MODEL, treeModel.getPrefix() + TreeModel.TREE_MODEL);
        writer.writeCloseTag(TMRCAStatisticParser.TMRCA_STATISTIC);
    }

    public void writeTMRCAStatisticReferences(XMLWriter writer) {
        for (Taxa taxa : options.taxonSets) {
            // make tmrca(tree.name) eay to read in log for Tracer
            PartitionTreeModel treeModel = options.taxonSetsTreeModel.get(taxa);
            writer.writeIDref(TMRCAStatisticParser.TMRCA_STATISTIC, "tmrca(" + treeModel.getPrefix() + taxa.getId() + ")");

        }
        for (Taxa taxa : options.taxonSets) {
            // make tmrca(tree.name) eay to read in log for Tracer
            PartitionTreeModel treeModel = options.taxonSetsTreeModel.get(taxa);
            if (treeModel.hasTipCalibrations()) {
                writer.writeIDref(TMRCAStatisticParser.TMRCA_STATISTIC, "age(" + treeModel.getPrefix() + taxa.getId() + ")");
            }
        }
    }

}