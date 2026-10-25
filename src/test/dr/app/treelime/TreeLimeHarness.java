/*
 * TreeLimeHarness.java
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

package test.dr.app.treelime;

import dr.app.treelime.components.ComponentFactory;
import dr.app.treelime.components.ancestralstates.AncestralStatesComponentFactory;
import dr.app.treelime.components.continuous.ContinuousComponentFactory;
import dr.app.treelime.components.discrete.DiscreteTraitsComponentFactory;
import dr.app.treelime.components.dollo.DolloComponentFactory;
import dr.app.treelime.components.hpm.HierarchicalModelComponentFactory;
import dr.app.treelime.components.linkedparameters.LinkedParameterComponentFactory;
import dr.app.treelime.components.marginalLikelihoodEstimation.MarginalLikelihoodEstimationComponentFactory;
import dr.app.treelime.components.sequenceerror.SequenceErrorModelComponentFactory;
import dr.app.treelime.components.tipdatesampling.TipDateSamplingComponentFactory;
import dr.app.beast.BeastParser;
import dr.app.beast.BeastVersion;
import dr.app.treelime.generator.BeastGenerator;
import dr.app.treelime.options.AbstractPartitionData;
import dr.app.treelime.options.PartitionClockModel;
import dr.app.treelime.options.PartitionTreeModel;
import dr.app.treelime.options.TreeLimeOptions;
import dr.app.treelime.util.TreeLimeImporter;
import dr.evolution.util.Date;
import dr.evolution.util.Taxon;
import dr.evolution.util.Units;
import dr.xml.XMLObject;
import junit.framework.Assert;

import java.io.File;
import java.io.StringReader;
import java.nio.file.Files;
import java.util.HashMap;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Drives TreeLime's generator without the GUI: import data, let a test set the options, write the XML, and optionally
 * check its ids and parse it with BEAST.
 */
class TreeLimeHarness {

    static final String TWO_PARTITIONS = "examples/Data/Primates_two_partition.nex";

    static final ComponentFactory[] COMPONENTS = {
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
    static TreeLimeOptions importData(String path) throws Exception {
        File data = new File(path);
        if (!data.exists()) {
            System.err.println("Skipping: " + path + " not found from " + new File(".").getAbsolutePath());
            return null;
        }
        TreeLimeOptions options = new TreeLimeOptions(COMPONENTS);
        new TreeLimeImporter(null, options).importFromFile(data);
        return options;
    }

    // dates from the two-digit year that ends each name (D4Brazi82 is 1982), as the Tip Dates panel would set them
    static void setTipDatesFromNames(TreeLimeOptions options) {
        Date mostRecent = null;
        for (int i = 0; i < options.taxonList.getTaxonCount(); i++) {
            Taxon taxon = options.taxonList.getTaxon(i);
            String name = taxon.getId();
            Date date = new Date(1900 + Integer.parseInt(name.substring(name.length() - 2)), Units.Type.YEARS, false);
            taxon.setDate(date);
            if (mostRecent == null || date.after(mostRecent)) {
                mostRecent = date;
            }
        }
        options.useTipDates = true;
        options.maximumTipHeight = 0.0;
        for (int i = 0; i < options.taxonList.getTaxonCount(); i++) {
            Taxon taxon = options.taxonList.getTaxon(i);
            double height = mostRecent.getTimeValue() - taxon.getDate().getTimeValue();
            taxon.setAttribute("_height", height);
            options.maximumTipHeight = Math.max(options.maximumTipHeight, height);
        }
    }

    // as the Data panel's "unlink clock models" does
    static void unlinkClockModels(TreeLimeOptions options) {
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
    static void unlinkTreeModels(TreeLimeOptions options) {
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
     * Writes the XML, checks that every idref resolves and no id repeats, and parses it with a one-iteration
     * maximizer; returns the XML text.
     */
    static String generateAndParse(TreeLimeOptions options) throws Exception {
        options.maximizationOptions.nIterations = 1;
        options.maximizationOptions.printToScreen = false;
        String xml = generate(options);
        assertIdsResolve(xml);

        // the reports would write files: drop them
        String parsed = xml.substring(0, xml.indexOf("<!-- The optimum")) + "</beast>\n";
        BeastParser parser = new BeastParser(new String[]{}, null, false, false, false, new BeastVersion());
        parser.parse(new StringReader(parsed), XMLObject.class);
        return xml;
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

    static String generate(TreeLimeOptions options) throws Exception {
        File xml = File.createTempFile("treeLime", ".xml");
        xml.deleteOnExit();
        options.fileNameStem = xml.getName().replace(".xml", "");
        options.selectParameters();
        options.selectOperators();
        new BeastGenerator(options, COMPONENTS).generateXML(xml);
        return new String(Files.readAllBytes(xml.toPath()));
    }
}
