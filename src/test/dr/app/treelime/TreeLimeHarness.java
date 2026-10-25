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
import dr.app.treelime.generator.BeastGenerator;
import dr.app.treelime.options.TreeLimeOptions;
import dr.app.treelime.util.TreeLimeImporter;

import java.io.File;
import java.nio.file.Files;

/**
 * Drives TreeLime's generator without the GUI: import data, let a test set the options, write the XML.
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
