/*
 * MixedEffectsClockGeneratorTest.java
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
import dr.app.beauti.types.ClockType;
import dr.app.beauti.util.BEAUTiImporter;
import dr.evolution.util.Taxa;
import dr.xml.XMLObject;
import junit.framework.TestCase;

import java.io.File;
import java.io.FileReader;
import java.nio.file.Files;
import java.util.HashMap;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * The mixed-effects clock XML that BEAUti writes must parse (plan 3, E7): one fixed effect and one category per
 * monophyletic taxon set, consecutively numbered, and every id prefixed so that several clock models can coexist.
 */
public class MixedEffectsClockGeneratorTest extends TestCase {

    private static final String DATA = "examples/Data/Primates_two_partition.nex";

    public MixedEffectsClockGeneratorTest(String name) {
        super(name);
    }

    public void testOneClockWithNonMonophyleticTaxonSet() throws Exception {
        checkGeneratedXmlParses(false);
    }

    public void testTwoUnlinkedClocks() throws Exception {
        checkGeneratedXmlParses(true);
    }

    private void checkGeneratedXmlParses(boolean unlinkClocks) throws Exception {
        File data = new File(DATA);
        if (!data.exists()) {
            System.err.println("Skipping: " + DATA + " not found from " + new File(".").getAbsolutePath());
            return;
        }

        ComponentFactory[] components = {
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
        BeautiOptions options = new BeautiOptions(components);
        new BEAUTiImporter(null, options).importFromFile(data);

        File xml = File.createTempFile("mixedEffectsClock", ".xml");
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

        // a non-monophyletic set first, so the numbering of the monophyletic ones is exercised
        int[][] sets = {{0, 7}, {1, 2}, {3, 4}};
        boolean[] mono = {false, true, true};
        for (int i = 0; i < sets.length; i++) {
            Taxa taxa = new Taxa("set" + (i + 1));
            for (int index : sets[i]) {
                taxa.addTaxon(options.taxonList.getTaxon(index));
            }
            options.taxonSets.add(taxa);
            options.taxonSetsMono.put(taxa, mono[i]);
            options.taxonSetsIncludeStem.put(taxa, i == 2);
            options.taxonSetsHeights.put(taxa, 0.0);
            options.taxonSetsTreeModel.put(taxa, options.getPartitionTreeModels().get(0));
        }

        if (unlinkClocks) {
            for (AbstractPartitionData partition : options.dataPartitions) {
                PartitionClockModel clockModel = partition.getPartitionClockModel();
                if (!clockModel.getName().equals(partition.getName())) {
                    clockModel = new PartitionClockModel(options, partition.getName(), clockModel);
                    partition.setPartitionClockModel(clockModel);
                }
                clockModel.setPartitionTreeModel(partition.getPartitionTreeModel());
            }
            assertEquals(2, options.getPartitionClockModels().size());
        }
        for (PartitionClockModel clock : options.getPartitionClockModels()) {
            clock.setClockType(ClockType.MIXED_EFFECTS_CLOCK);
        }
        options.selectParameters();
        options.selectOperators();
        new BeastGenerator(options, components).generateXML(xml);

        String text = new String(Files.readAllBytes(xml.toPath()));
        assertNoDanglingIdrefs(text);
        assertTrue(text.contains("usePreOrder=\"true\""));

        // parsing builds the fixed-effects model, which checks the coefficient dimension
        BeastParser parser = new BeastParser(new String[]{}, null, false, false, false, new BeastVersion());
        parser.parse(new FileReader(xml), XMLObject.class);
    }

    private static void assertNoDanglingIdrefs(String xml) {
        Map<String, Integer> ids = new HashMap<>();
        Matcher id = Pattern.compile("\\sid=\"([^\"]+)\"").matcher(xml);
        while (id.find()) {
            ids.merge(id.group(1), 1, Integer::sum);
        }
        for (Map.Entry<String, Integer> entry : ids.entrySet()) {
            assertEquals("duplicate id " + entry.getKey(), 1, (int) entry.getValue());
        }
        Matcher idref = Pattern.compile("idref=\"([^\"]+)\"").matcher(xml);
        while (idref.find()) {
            assertTrue("dangling idref " + idref.group(1), ids.containsKey(idref.group(1)));
        }
    }
}
