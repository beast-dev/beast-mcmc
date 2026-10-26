/*
 * UnsupportedGradientRefusalTest.java
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

package test.dr.evomodel.tree;

import dr.app.beast.BeastParser;
import dr.app.beast.BeastVersion;
import dr.evomodelxml.tree.SampledTipHeights;
import dr.inference.model.Parameter;
import junit.framework.TestCase;

import java.io.StringReader;

/**
 * Gradients with respect to sampled tip dates (plan 3, F2) and for the old birth-death serial-sampling model (G7)
 * are not implemented, so building one must fail with a clear message.
 */
public class UnsupportedGradientRefusalTest extends TestCase {

    private static final String TREE =
            "  <taxa id=\"taxa\"><taxon id=\"a\"/><taxon id=\"b\"/><taxon id=\"c\"/><taxon id=\"d\"/></taxa>\n" +
            "  <newick id=\"startingTree\">(((a:1,b:1):1,c:2):1,d:3);</newick>\n" +
            "  <treeModel id=\"treeModel\">\n" +
            "    <tree idref=\"startingTree\"/>\n" +
            "    <rootHeight><parameter id=\"treeModel.rootHeight\"/></rootHeight>\n" +
            "    <nodeHeights internalNodes=\"true\"><parameter id=\"treeModel.internalNodeHeights\"/></nodeHeights>\n" +
            "    <nodeHeights internalNodes=\"true\" rootNode=\"true\">" +
            "<parameter id=\"treeModel.allInternalNodeHeights\"/></nodeHeights>\n" +
            "    LEAF\n" +
            "  </treeModel>\n";

    public UnsupportedGradientRefusalTest(String name) {
        super(name);
    }

    public void testGradientWrtSampledTipDate() {
        String xml = "<beast>\n" +
                TREE.replace("LEAF", "<leafHeight taxon=\"a\"><parameter id=\"age(a)\"/></leafHeight>") +
                "  <uniformPrior id=\"agePrior\" lower=\"0.0\" upper=\"0.5\"><parameter idref=\"age(a)\"/></uniformPrior>\n" +
                "  <gradient id=\"gradient\"><uniformPrior idref=\"agePrior\"/><parameter idref=\"age(a)\"/></gradient>\n" +
                "</beast>\n";
        assertRefused(xml, "sampled tip dates");
    }

    public void testGradientWithoutSampledTipDateStillParses() throws Exception {
        String xml = "<beast>\n" + TREE.replace("LEAF", "") +
                "  <parameter id=\"x\" value=\"0.2\" lower=\"0.0\"/>\n" +
                "  <uniformPrior id=\"xPrior\" lower=\"0.0\" upper=\"0.5\"><parameter idref=\"x\"/></uniformPrior>\n" +
                "  <gradient id=\"gradient\"><uniformPrior idref=\"xPrior\"/><parameter idref=\"x\"/></gradient>\n" +
                "</beast>\n";
        parser().parse(new StringReader(xml), null);
    }

    // the trees are found among the parsed objects; a leaf height counts alone or inside a compound parameter
    public void testSampledTipHeightLookup() throws Exception {
        String xml = "<beast>\n" +
                TREE.replace("LEAF", "<leafHeight taxon=\"a\"><parameter id=\"age(a)\"/></leafHeight>") +
                "  <parameter id=\"x\" value=\"0.2\"/>\n" +
                "  <compoundParameter id=\"both\"><parameter idref=\"x\"/><parameter idref=\"age(a)\"/></compoundParameter>\n" +
                "</beast>\n";
        java.util.Map<String, dr.xml.XMLObject> store = parser().parse(new StringReader(xml), false);
        assertTrue(SampledTipHeights.containsSampledTipHeight(parameter(store, "age(a)"), store));
        assertTrue(SampledTipHeights.containsSampledTipHeight(parameter(store, "both"), store));
        assertFalse(SampledTipHeights.containsSampledTipHeight(parameter(store, "x"), store));
        assertFalse(SampledTipHeights.containsSampledTipHeight(parameter(store, "treeModel.rootHeight"), store));
        assertFalse(SampledTipHeights.containsSampledTipHeight(parameter(store, "treeModel.allInternalNodeHeights"),
                store));
    }

    private static Parameter parameter(java.util.Map<String, dr.xml.XMLObject> store, String id) {
        return (Parameter) store.get(id).getNativeObject();
    }

    public void testBirthDeathSerialSampling() {
        String xml = "<beast>\n" + TREE.replace("LEAF", "") +
                "  <birthDeathSerialSampling id=\"bdss\" units=\"years\">\n" +
                "    <birthRate><parameter id=\"bdss.birthRate\" value=\"2.0\"/></birthRate>\n" +
                "    <relativeDeathRate><parameter id=\"bdss.relativeDeathRate\" value=\"0.5\"/></relativeDeathRate>\n" +
                "    <psi><parameter id=\"bdss.psi\" value=\"0.5\"/></psi>\n" +
                "    <origin><parameter id=\"bdss.origin\" value=\"5.0\"/></origin>\n" +
                "  </birthDeathSerialSampling>\n" +
                "  <speciationLikelihood id=\"speciation\">\n" +
                "    <model><birthDeathSerialSampling idref=\"bdss\"/></model>\n" +
                "    <speciesTree><treeModel idref=\"treeModel\"/></speciesTree>\n" +
                "  </speciationLikelihood>\n" +
                "  <speciationLikelihoodGradient id=\"gradient\" wrtParameter=\"nodeHeight\">" +
                "<speciationLikelihood idref=\"speciation\"/><treeModel idref=\"treeModel\"/>" +
                "</speciationLikelihoodGradient>\n" +
                "</beast>\n";
        assertRefused(xml, "BDSS");
    }

    private static BeastParser parser() {
        return new BeastParser(new String[]{}, null, false, false, false, new BeastVersion());
    }

    private static void assertRefused(String xml, String expectedMessage) {
        try {
            parser().parse(new StringReader(xml), null);
            fail("expected the gradient to be refused");
        } catch (Exception e) {
            String message = e.toString();
            Throwable cause = e.getCause();
            while (cause != null) {
                message += " / " + cause;
                cause = cause.getCause();
            }
            assertTrue(message, message.contains(expectedMessage));
        }
    }
}
