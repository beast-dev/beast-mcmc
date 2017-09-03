/*
 * TreeModelGenerator.java
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
import dr.app.beauti.options.AbstractPartitionData;
import dr.app.beauti.options.BeautiOptions;
import dr.app.beauti.options.PartitionTreeModel;
import dr.app.beauti.util.XMLWriter;
import dr.evolution.datatype.DataType;
import dr.evomodel.tree.TMRCAStatistic;
import dr.evomodel.tree.TreeModel;
import dr.evomodelxml.coalescent.OldCoalescentSimulatorParser;
import dr.evomodelxml.tree.MicrosatelliteSamplerTreeModelParser;
import dr.evomodelxml.tree.TMRCAStatisticParser;
import dr.evomodelxml.tree.TreeModelParser;
import dr.evoxml.MicrosatellitePatternParser;
import dr.evoxml.UPGMATreeParser;
import dr.inference.model.ParameterParser;
import dr.util.Attribute;
import dr.xml.XMLParser;

/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @author Walter Xie
 */
public class TreeModelGenerator extends Generator {

    public TreeModelGenerator(BeautiOptions options, ComponentFactory[] components) {
        super(options, components);
    }

    /**
     * Write tree model XML block.
     *
     * @param model
     * @param writer the writer
     */
    void writeTreeModel(PartitionTreeModel model, XMLWriter writer) {

        String prefix = model.getPrefix();

        final String treeModelName = prefix + TreeModel.TREE_MODEL; // treemodel.treeModel or treeModel

        writer.writeComment("Generate a tree model");
        writer.writeTag(TreeModel.TREE_MODEL, new Attribute.Default<String>(XMLParser.ID, treeModelName), false);

        final String STARTING_TREE = InitialTreeGenerator.STARTING_TREE;

        switch (model.getStartingTreeType()) {
            case USER:
                writer.writeIDref("tree", prefix + STARTING_TREE);
                break;
            case UPGMA:
                writer.writeIDref(UPGMATreeParser.UPGMA_TREE, prefix + STARTING_TREE);
                break;
            case RANDOM:
                writer.writeIDref(OldCoalescentSimulatorParser.COALESCENT_TREE, prefix + STARTING_TREE);
                break;
            default:
                throw new IllegalArgumentException("Unknown StartingTreeType");
        }

        writer.writeOpenTag(TreeModelParser.ROOT_HEIGHT);
        writer.writeTag(ParameterParser.PARAMETER,
                new Attribute.Default<String>(XMLParser.ID, treeModelName + "." + OldCoalescentSimulatorParser.ROOT_HEIGHT), true);
        writer.writeCloseTag(TreeModelParser.ROOT_HEIGHT);


        writer.writeOpenTag(TreeModelParser.NODE_HEIGHTS, new Attribute.Default<String>(TreeModelParser.INTERNAL_NODES, "true"));
        writer.writeTag(ParameterParser.PARAMETER,
                new Attribute.Default<String>(XMLParser.ID, treeModelName + "." + "internalNodeHeights"), true);
        writer.writeCloseTag(TreeModelParser.NODE_HEIGHTS);

        writer.writeOpenTag(TreeModelParser.NODE_HEIGHTS,
                new Attribute[]{
                        new Attribute.Default<String>(TreeModelParser.INTERNAL_NODES, "true"),
                        new Attribute.Default<String>(TreeModelParser.ROOT_NODE, "true")
                });
        writer.writeTag(ParameterParser.PARAMETER,
                new Attribute.Default<String>(XMLParser.ID, treeModelName + "." + "allInternalNodeHeights"), true);
        writer.writeCloseTag(TreeModelParser.NODE_HEIGHTS);

//        int randomLocalClockCount = 0;
//        int autocorrelatedClockCount = 0;
//        for (PartitionData pd : model.getDataPartitions()) { // only the PDs linked to this tree model
//        	PartitionClockModel clockModel = pd.getPartitionClockModel();
//        	switch (clockModel.getClockType()) {
//	        	case AUTOCORRELATED_LOGNORMAL: autocorrelatedClockCount += 1; break;
//	        	case RANDOM_LOCAL_CLOCK: randomLocalClockCount += 1; break;
//        	}
//        }
//
//        if (autocorrelatedClockCount > 1 || randomLocalClockCount > 1 || autocorrelatedClockCount + randomLocalClockCount > 1) {
//        	//FAIL
//            throw new IllegalArgumentException("clock model/tree model combination not implemented by BEAST yet!");
//        }
        // move to validateClockTreeModelCombination(PartitionTreeModel model)

//    	if (autocorrelatedClockCount == 1) {
//        if (count[0] == 1) {
//                writer.writeOpenTag(TreeModelParser.NODE_RATES,
//                        new Attribute[]{
//                                new Attribute.Default<String>(TreeModelParser.ROOT_NODE, "false"),
//                                new Attribute.Default<String>(TreeModelParser.INTERNAL_NODES, "true"),
//                                new Attribute.Default<String>(TreeModelParser.LEAF_NODES, "true")
//                        });
//                writer.writeTag(ParameterParser.PARAMETER,
//                        new Attribute.Default<String>(XMLParser.ID, treeModelName + "." + TreeModelParser.NODE_RATES), true);
//                writer.writeCloseTag(TreeModelParser.NODE_RATES);
//
//                writer.writeOpenTag(TreeModelParser.NODE_RATES,
//                        new Attribute[]{
//                                new Attribute.Default<String>(TreeModelParser.ROOT_NODE, "true"),
//                                new Attribute.Default<String>(TreeModelParser.INTERNAL_NODES, "false"),
//                                new Attribute.Default<String>(TreeModelParser.LEAF_NODES, "false")
//                        });
//                writer.writeTag(ParameterParser.PARAMETER,
//                        new Attribute.Default<String>(XMLParser.ID,
//                                treeModelName + "." + RateEvolutionLikelihood.ROOTRATE), true);
//                writer.writeCloseTag(TreeModelParser.NODE_RATES);
////    	} else if (randomLocalClockCount == 1 ) {
//        } else

        //+++++++++++++ removed because random local clock XML is changed ++++++++++++++++
//        int[] count = validateClockTreeModelCombination(model);
//        if (count[1] == 1) {
//            writer.writeOpenTag(TreeModelParser.NODE_RATES,
//                    new Attribute[]{
//                            new Attribute.Default<String>(TreeModelParser.ROOT_NODE, "false"),
//                            new Attribute.Default<String>(TreeModelParser.INTERNAL_NODES, "true"),
//                            new Attribute.Default<String>(TreeModelParser.LEAF_NODES, "true")
//                    });
//            writer.writeTag(ParameterParser.PARAMETER,
//                    new Attribute.Default<String>(XMLParser.ID, modelPrefix + ClockType.LOCAL_CLOCK + ".relativeRates"), true);
//            writer.writeCloseTag(TreeModelParser.NODE_RATES);
//
//            writer.writeOpenTag(TreeModelParser.NODE_TRAITS,
//                    new Attribute[]{
//                            new Attribute.Default<String>(TreeModelParser.ROOT_NODE, "false"),
//                            new Attribute.Default<String>(TreeModelParser.INTERNAL_NODES, "true"),
//                            new Attribute.Default<String>(TreeModelParser.LEAF_NODES, "true")
//                    });
//            writer.writeTag(ParameterParser.PARAMETER,
//                    new Attribute.Default<String>(XMLParser.ID, modelPrefix + ClockType.LOCAL_CLOCK + ".changes"), true);
//            writer.writeCloseTag(TreeModelParser.NODE_TRAITS);
//        }

        generateInsertionPoint(ComponentGenerator.InsertionPoint.IN_TREE_MODEL, model, writer);

        writer.writeCloseTag(TreeModel.TREE_MODEL);

//        if (autocorrelatedClockCount == 1) {
//        if (count[0] == 1) {
//            writer.writeText("");
//            writer.writeOpenTag(CompoundParameter.COMPOUND_PARAMETER,
//                    new Attribute[]{new Attribute.Default<String>(XMLParser.ID, treeModelName + "." + "allRates")});
//            writer.writeIDref(ParameterParser.PARAMETER, treeModelName + "." + TreeModelParser.NODE_RATES);
//            writer.writeIDref(ParameterParser.PARAMETER, treeModelName + "." + RateEvolutionLikelihood.ROOTRATE);
//            writer.writeCloseTag(CompoundParameter.COMPOUND_PARAMETER);
//        }

        if (model.hasTipCalibrations()) {
            writer.writeComment("Statistic for time of most recent common ancestor of tree");
            writer.writeTag(TMRCAStatisticParser.TMRCA_STATISTIC,
                    new Attribute[]{
                            new Attribute.Default<String>(XMLParser.ID, prefix + "age(root)"),
                            new Attribute.Default<String>(TMRCAStatisticParser.ABSOLUTE, "true")
                    }, false);
            writer.writeIDref(TreeModel.TREE_MODEL, treeModelName);
            writer.writeCloseTag(TMRCAStatisticParser.TMRCA_STATISTIC);
        }

        if (model.getDataType().getType() == DataType.MICRO_SAT) {
            for (AbstractPartitionData partitionData : options.getDataPartitions(model)) {
                writer.writeComment("Generate a microsatellite tree model");
                writer.writeTag(MicrosatelliteSamplerTreeModelParser.TREE_MICROSATELLITE_SAMPLER_MODEL,
                        new Attribute.Default<String>(XMLParser.ID, partitionData.getName() + "." +
                                MicrosatelliteSamplerTreeModelParser.TREE_MICROSATELLITE_SAMPLER_MODEL), false);

                writer.writeOpenTag(MicrosatelliteSamplerTreeModelParser.TREE);
                writer.writeIDref(TreeModel.TREE_MODEL, treeModelName);
                writer.writeCloseTag(MicrosatelliteSamplerTreeModelParser.TREE);

                writer.writeOpenTag(MicrosatelliteSamplerTreeModelParser.INTERNAL_VALUES);
                writer.writeTag(ParameterParser.PARAMETER, new Attribute[]{
                        new Attribute.Default<String>(XMLParser.ID, partitionData.getName() + "." +
                                MicrosatelliteSamplerTreeModelParser.TREE_MICROSATELLITE_SAMPLER_MODEL + ".internalNodesParameter"),
                        new Attribute.Default<Integer>(ParameterParser.DIMENSION, model.getDimension())}, true);
                writer.writeCloseTag(MicrosatelliteSamplerTreeModelParser.INTERNAL_VALUES);

                writer.writeOpenTag(MicrosatelliteSamplerTreeModelParser.EXTERNAL_VALUES);

                writer.writeIDref(MicrosatellitePatternParser.MICROSATPATTERN, partitionData.getName());

                writer.writeCloseTag(MicrosatelliteSamplerTreeModelParser.EXTERNAL_VALUES);

                writer.writeCloseTag(MicrosatelliteSamplerTreeModelParser.TREE_MICROSATELLITE_SAMPLER_MODEL);
            }
        }
    }
}
