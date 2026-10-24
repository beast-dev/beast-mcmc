/*
 * BeautiGeneratorHarness.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
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
 *
 */

package test.dr.app.beauti;

import dr.app.beast.BeastParser;
import dr.app.beast.BeastVersion;
import dr.app.beauti.components.ComponentFactory;
import dr.app.beauti.components.ancestralstates.AncestralStatesComponentFactory;
import dr.app.beauti.components.continuous.ContinuousComponentFactory;
import dr.app.beauti.components.discrete.DiscreteTraitsComponentFactory;
import dr.app.beauti.components.dollo.DolloComponentFactory;
import dr.app.beauti.components.hpm.HierarchicalModelComponentFactory;
import dr.app.beauti.components.linkedparameters.LinkedParameterComponentFactory;
import dr.app.beauti.components.marginalLikelihoodEstimation.MarginalLikelihoodEstimationComponentFactory;
import dr.app.beauti.components.sequenceerror.SequenceErrorModelComponentFactory;
import dr.app.beauti.components.tipdatesampling.TipDateSamplingComponentFactory;
import dr.app.beauti.generator.BeastGenerator;
import dr.app.beauti.options.AbstractPartitionData;
import dr.app.beauti.options.BeautiOptions;
import dr.app.beauti.options.PartitionClockModel;
import dr.app.beauti.options.PartitionTreeModel;
import dr.app.beauti.util.BEAUTiImporter;
import dr.xml.XMLObject;
import junit.framework.Assert;

import java.io.File;
import java.io.FileReader;
import java.nio.file.Files;
import java.util.HashMap;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Drives BEAUti's generator without the GUI: import data, let a test set the model options, write the XML, check
 * that every idref resolves and no id repeats, and parse it with BEAST.
 */
class BeautiGeneratorHarness {

    static final String TWO_PARTITIONS = "examples/Data/Primates_two_partition.nex";

    private static final ComponentFactory[] COMPONENTS = {
            AncestralStatesComponentFactory.INSTANCE,
            ContinuousComponentFactory.INSTANCE,
            DiscreteTraitsComponentFactory.INSTANCE,
            DolloComponentFactory.INSTANCE,
            LinkedParameterComponentFactory.INSTANCE,
            HierarchicalModelComponentFactory.INSTANCE,
            MarginalLikelihoodEstimationComponentFactory.INSTANCE,
            SequenceErrorModelComponentFactory.INSTANCE,
            TipDateSamplingComponentFactory.INSTANCE
    };

    /**
     * @return the options with the data imported, or null (with a message) if the data file is not found
     */
    static BeautiOptions importData(String path) throws Exception {
        File data = new File(path);
        if (!data.exists()) {
            System.err.println("Skipping: " + path + " not found from " + new File(".").getAbsolutePath());
            return null;
        }
        BeautiOptions options = new BeautiOptions(COMPONENTS);
        new BEAUTiImporter(null, options).importFromFile(data);
        return options;
    }

    // as the Data panel's "unlink clock models" does
    static void unlinkClockModels(BeautiOptions options) {
        for (AbstractPartitionData partition : options.dataPartitions) {
            PartitionClockModel clockModel = partition.getPartitionClockModel();
            if (!clockModel.getName().equals(partition.getName())) {
                clockModel = new PartitionClockModel(options, partition.getName(), clockModel);
                partition.setPartitionClockModel(clockModel);
            }
            clockModel.setPartitionTreeModel(partition.getPartitionTreeModel());
        }
    }

    // as the Data panel's "unlink trees" does (which also unlinks the clocks and links the tree priors)
    static void unlinkTreeModels(BeautiOptions options) {
        for (AbstractPartitionData partition : options.dataPartitions) {
            PartitionTreeModel treeModel = partition.getPartitionTreeModel();
            if (!treeModel.getName().equals(partition.getName())) {
                treeModel = new PartitionTreeModel(options, partition.getName(), treeModel);
                partition.setPartitionTreeModel(treeModel);
            }
            PartitionClockModel clockModel = partition.getPartitionClockModel();
            if (!clockModel.getName().equals(partition.getName())) {
                clockModel = new PartitionClockModel(options, partition.getName(), clockModel);
                partition.setPartitionClockModel(clockModel);
            }
            clockModel.setPartitionTreeModel(treeModel);
        }
        options.linkTreePriors(options.getPartitionTreeModels().get(0).getPartitionTreePrior());
    }

    /**
     * Writes the XML, checks its ids and parses it; returns the XML text.
     */
    static String generateAndParse(BeautiOptions options) throws Exception {
        File xml = File.createTempFile("beautiGenerator", ".xml");
        xml.deleteOnExit();
        options.fileNameStem = xml.getName().replace(".xml", "");
        options.logFileName = options.fileNameStem + ".log";
        options.chainLength = 10;
        options.logEvery = 10;
        options.echoEvery = 10;
        options.treeFileName.clear();
        for (PartitionTreeModel tree : options.getPartitionTreeModels()) {
            options.treeFileName.add(options.fileNameStem + "." + tree.getPrefix() + "trees");
        }
        options.selectParameters();
        options.selectOperators();
        new BeastGenerator(options, COMPONENTS).generateXML(xml);

        String text = new String(Files.readAllBytes(xml.toPath()));
        assertIdsResolve(text);

        BeastParser parser = new BeastParser(new String[]{}, null, false, false, false, new BeastVersion());
        parser.parse(new FileReader(xml), XMLObject.class);
        return text;
    }

    static void assertIdsResolve(String xml) {
        Map<String, Integer> ids = new HashMap<>();
        Matcher id = Pattern.compile("\\sid=\"([^\"]+)\"").matcher(xml);
        while (id.find()) {
            ids.merge(id.group(1), 1, Integer::sum);
        }
        for (Map.Entry<String, Integer> entry : ids.entrySet()) {
            Assert.assertEquals("duplicate id " + entry.getKey(), 1, (int) entry.getValue());
        }
        Matcher idref = Pattern.compile("idref=\"([^\"]+)\"").matcher(xml);
        while (idref.find()) {
            Assert.assertTrue("dangling idref " + idref.group(1), ids.containsKey(idref.group(1)));
        }
    }
}
