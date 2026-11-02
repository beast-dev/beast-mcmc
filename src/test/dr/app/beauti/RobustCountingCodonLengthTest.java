/*
 * RobustCountingCodonLengthTest.java
 *
 * Copyright © 2002-2026 the BEAST Development Team
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
 */

package test.dr.app.beauti;

import dr.app.beauti.components.ancestralstates.AncestralStatesComponentOptions;
import dr.app.beauti.generator.BeastGenerator;
import dr.app.beauti.generator.Generator;
import dr.app.beauti.options.AbstractPartitionData;
import dr.app.beauti.options.BeautiOptions;
import dr.app.beauti.options.PartitionSubstitutionModel;
import dr.app.beauti.util.BeautiImporter;
import org.junit.Test;

import java.io.File;
import java.io.PrintWriter;
import java.lang.reflect.Field;
import java.lang.reflect.Method;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;

import static org.junit.Assert.*;

/**
 * Renaissance (robust) counting compares the three codon positions site by site, so each must hold the same
 * number of sites. Taking every third site of a partition gives three equal counts only when the partition
 * is a whole number of codons long; otherwise BEAST refuses the file with "Codon-partitioned robust counting
 * requires all partitions to have the same length".
 *
 * BEAUti offered the option without checking this, so the problem only showed up when the run was started.
 */
public class RobustCountingCodonLengthTest {

    /**
     * The two-partition primates example, with the two charsets trimmed to the given length so that the
     * partitions are, or are not, a whole number of codons long.
     */
    private File dataWithPartitionsOfLength(int length) throws Exception {
        File source = new File(repositoryRoot(), "examples/Data/Primates_two_partition.nex");
        String nexus = new String(Files.readAllBytes(source.toPath()), StandardCharsets.UTF_8);

        nexus = nexus.replace("charset firsthalf = 1-449;", "charset firsthalf = 1-" + length + ";");
        nexus = nexus.replace("charset secondhalf = 450-898;",
                "charset secondhalf = " + (length + 1) + "-" + (2 * length) + ";");

        File file = File.createTempFile("primates", ".nex");
        file.deleteOnExit();
        PrintWriter writer = new PrintWriter(file, "UTF-8");
        writer.print(nexus);
        writer.close();
        return file;
    }

    private File repositoryRoot() {
        File dir = new File(System.getProperty("user.dir"));
        while (dir != null && !new File(dir, "examples/Data").isDirectory()) {
            dir = dir.getParentFile();
        }
        assertNotNull("could not find the repository root", dir);
        return dir;
    }

    /** the two-partition primates data with three codon positions per partition */
    private BeautiOptions options(int partitionLength) throws Exception {
        BeautiOptions options = BeautiOptions.newInstance(BeautiHeadlessHarness.COMPONENTS);
        new BeautiImporter(null, options).importFromFile(dataWithPartitionsOfLength(partitionLength));
        BeautiHeadlessHarness.unlinkSubstitutionModels(options);
        for (PartitionSubstitutionModel model : options.getPartitionSubstitutionModels()) {
            model.setCodonHeteroPattern("123");
        }
        return options;
    }

    private AncestralStatesComponentOptions ancestralStates(BeautiOptions options) {
        return (AncestralStatesComponentOptions)
                options.getComponentOptions(AncestralStatesComponentOptions.class);
    }

    /** sets the flag behind the setter's back, as loading a document saved by an older BEAUti would */
    private void forceRobustCountingOn(BeautiOptions options, AbstractPartitionData partition) throws Exception {
        Method getOptions = AncestralStatesComponentOptions.class
                .getDeclaredMethod("getOptions", AbstractPartitionData.class);
        getOptions.setAccessible(true);
        Object inner = getOptions.invoke(ancestralStates(options), partition);
        Field field = inner.getClass().getDeclaredField("dNdSRobustCounting");
        field.setAccessible(true);
        field.setBoolean(inner, true);
    }

    /** @return the message if generation was refused, or null if it succeeded */
    private String checkOptions(BeautiOptions options) {
        try {
            new BeastGenerator(options, BeautiHeadlessHarness.COMPONENTS).checkOptions();
            return null;
        } catch (Generator.GeneratorException e) {
            return e.getMessage();
        }
    }

    @Test
    public void aWholeNumberOfCodonsIsAccepted() throws Exception {
        BeautiOptions options = options(447);

        for (AbstractPartitionData partition : options.dataPartitions) {
            assertEquals(447, partition.getSiteCount());
            assertTrue("447 sites is 149 codons", ancestralStates(options).isWholeNumberOfCodons(partition));
            assertTrue(ancestralStates(options).dNdSRobustCountingAvailable(partition));
            assertNull(ancestralStates(options).dNdSRobustCountingUnavailableReason(partition));
            ancestralStates(options).setDNdSRobustCounting(partition, true);
            assertTrue(ancestralStates(options).dNdSRobustCounting(partition));
        }

        assertNull("generation should not be refused", checkOptions(options));
    }

    @Test
    public void aPartLastCodonMakesTheOptionUnavailable() throws Exception {
        BeautiOptions options = options(449);

        for (AbstractPartitionData partition : options.dataPartitions) {
            assertEquals(449, partition.getSiteCount());
            assertFalse(ancestralStates(options).isWholeNumberOfCodons(partition));
            assertFalse(ancestralStates(options).dNdSRobustCountingAvailable(partition));

            String reason = ancestralStates(options).dNdSRobustCountingUnavailableReason(partition);
            assertNotNull(reason);
            assertTrue(reason, reason.contains("449"));
            assertTrue(reason, reason.contains("codon"));
        }
    }

    @Test
    public void turningTheOptionOnForAPartLastCodonIsRefused() throws Exception {
        BeautiOptions options = options(449);
        AbstractPartitionData partition = options.dataPartitions.get(0);

        try {
            ancestralStates(options).setDNdSRobustCounting(partition, true);
            fail("should not be able to turn robust counting on for a partition that is not whole codons");
        } catch (IllegalArgumentException expected) {
            assertTrue(expected.getMessage(), expected.getMessage().contains("codon"));
        }

        assertFalse(ancestralStates(options).dNdSRobustCounting(partition));
    }

    @Test
    public void generationIsRefusedWhenTheOptionIsAlreadySet() throws Exception {
        // a document saved before this was checked would come back with the option set
        BeautiOptions options = options(449);
        for (AbstractPartitionData partition : options.dataPartitions) {
            forceRobustCountingOn(options, partition);
            assertTrue(ancestralStates(options).dNdSRobustCounting(partition));
        }

        String refused = checkOptions(options);
        assertNotNull("generation should have been refused", refused);
        assertTrue(refused, refused.contains("449"));
        assertTrue(refused, refused.contains(options.dataPartitions.get(0).getName()));
    }

    @Test
    public void theOtherRequirementsAreStillReported() throws Exception {
        BeautiOptions options = options(447);
        for (PartitionSubstitutionModel model : options.getPartitionSubstitutionModels()) {
            model.setCodonHeteroPattern(null); // no codon partitioning
        }

        AbstractPartitionData partition = options.dataPartitions.get(0);
        assertFalse(ancestralStates(options).dNdSRobustCountingAvailable(partition));
        String reason = ancestralStates(options).dNdSRobustCountingUnavailableReason(partition);
        assertNotNull(reason);
        assertTrue(reason, reason.contains("3-partition codon model"));
    }
}
