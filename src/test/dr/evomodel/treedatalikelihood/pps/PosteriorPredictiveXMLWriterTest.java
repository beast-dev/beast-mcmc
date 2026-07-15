/*
 * PosteriorPredictiveXMLWriterTest.java
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

package test.dr.evomodel.treedatalikelihood.pps;

import dr.evolution.alignment.SimpleAlignment;
import dr.evolution.datatype.DataType;
import dr.evolution.datatype.GeneralDataType;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.sequence.DelimitedSequence;
import dr.evolution.sequence.Sequence;
import dr.evolution.util.Taxa;
import dr.evolution.util.Taxon;
import dr.evomodel.treedatalikelihood.pps.PosteriorPredictiveXMLWriter;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

/**
 */
public class PosteriorPredictiveXMLWriterTest extends TestCase {

    private static final String NL = System.lineSeparator();

    private PosteriorPredictiveXMLWriter writer;

    public PosteriorPredictiveXMLWriterTest(String name) {
        super(name);
    }

    public void setUp() throws Exception {
        super.setUp();
        writer = new PosteriorPredictiveXMLWriter();
    }

    public void testWriteRootOpen() {
        assertEquals("<pps>" + NL, writer.writeRootOpen());
    }

    public void testWriteRootClose() {
        assertEquals("</pps>" + NL, writer.writeRootClose());
    }

    public void testWriteTaxaBlock() {
        Taxa taxa = new Taxa();
        taxa.addTaxon(new Taxon("A"));
        taxa.addTaxon(new Taxon("B"));

        String expected =
                "<taxa id=\"taxa\">" + NL +
                        "\t<taxon id=\"A\"/>" + NL +
                        "\t<taxon id=\"B\"/>" + NL +
                        "</taxa>" + NL;

        assertEquals(expected, writer.writeTaxaBlock(taxa, "taxa"));
    }

    public void testWriteGeneralDataTypeBlock() {
        GeneralDataType dataType = new GeneralDataType(new String[]{"location1", "location2", "location3"});
        dataType.setId("geography");

        String expected =
                "<generalDataType id=\"geography\">" + NL +
                        "\t<state code=\"location1\"/>" + NL +
                        "\t<state code=\"location2\"/>" + NL +
                        "\t<state code=\"location3\"/>" + NL +
                        "</generalDataType>" + NL;

        assertEquals(expected, writer.writeGeneralDataTypeBlock(dataType, "geography"));
    }

    public void testWriteAlignmentBlockForNucleotideData() {
        SimpleAlignment alignment = new SimpleAlignment();
        alignment.setDataType(Nucleotides.INSTANCE);
        addSequence(alignment, "A", "ACGT");
        addSequence(alignment, "B", "ACGA");

        String expected =
                "<alignment id=\"alignment_1\" dataType=\"nucleotide\" state=\"1000000\" parentId=\"patterns\">" + NL +
                        "\t<sequence>" + NL +
                        "\t\t<taxon idref=\"A\"/>" + NL +
                        "\t\tACGT" + NL +
                        "\t</sequence>" + NL +
                        "\t<sequence>" + NL +
                        "\t\t<taxon idref=\"B\"/>" + NL +
                        "\t\tACGA" + NL +
                        "\t</sequence>" + NL +
                        "</alignment>" + NL;

        assertEquals(expected, writer.writeAlignmentBlock(alignment, "alignment_1", 1000000L, "patterns"));
    }

    public void testWriteAlignmentBlockForGeneralData() {
        GeneralDataType dataType = new GeneralDataType(new String[]{"location1", "location2", "location3"});
        dataType.setId("geography");

        SimpleAlignment alignment = new SimpleAlignment();
        alignment.setDataType(dataType);
        addDelimitedSequence(alignment, dataType, "SimSeq1", "location2");
        addDelimitedSequence(alignment, dataType, "SimSeq2", "location3");

        String expected =
                "<alignment id=\"alignment_1\" dataType=\"geography\" state=\"0\" parentId=\"geoPatterns\">" + NL +
                        "\t<sequence>" + NL +
                        "\t\t<taxon idref=\"SimSeq1\"/>" + NL +
                        "\t\tlocation2" + NL +
                        "\t</sequence>" + NL +
                        "\t<sequence>" + NL +
                        "\t\t<taxon idref=\"SimSeq2\"/>" + NL +
                        "\t\tlocation3" + NL +
                        "\t</sequence>" + NL +
                        "</alignment>" + NL;

        assertEquals(expected, writer.writeAlignmentBlock(alignment, "alignment_1", 0L, "geoPatterns"));
    }

    private static void addSequence(SimpleAlignment alignment, String taxonId, String sequenceString) {
        Sequence sequence = new Sequence(sequenceString);
        sequence.setTaxon(new Taxon(taxonId));
        sequence.setDataType(alignment.getDataType());
        alignment.addSequence(sequence);
    }

    // GeneralDataType with multi-character state codes is delimited, and BeagleSequenceSimulator
    // represents delimited output with a case-preserving DelimitedSequence (see Utils.intArray2Sequence)
    // rather than the case-folding base Sequence.
    private static void addDelimitedSequence(SimpleAlignment alignment, DataType dataType,
                                              String taxonId, String sequenceString) {
        Sequence sequence = new DelimitedSequence(new Taxon(taxonId), sequenceString, dataType);
        alignment.addSequence(sequence);
    }

    public static Test suite() {
        return new TestSuite(PosteriorPredictiveXMLWriterTest.class);
    }
}
