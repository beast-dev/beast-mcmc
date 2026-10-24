/*
 * TMRCAPriorGradientRefusalTest.java
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
import junit.framework.TestCase;

import java.io.StringReader;

/**
 * Gradients with tMRCA calibration priors are not implemented (plan 3, F4), so building one must fail loudly instead
 * of silently leaving the calibration term out.
 */
public class TMRCAPriorGradientRefusalTest extends TestCase {

    private static final String TREE =
            "  <taxa id=\"taxa\"><taxon id=\"a\"/><taxon id=\"b\"/><taxon id=\"c\"/><taxon id=\"d\"/></taxa>\n" +
            "  <taxa id=\"clade\"><taxon idref=\"a\"/><taxon idref=\"b\"/></taxa>\n" +
            "  <newick id=\"startingTree\">(((a:1,b:1):1,c:2):1,d:3);</newick>\n" +
            "  <treeModel id=\"treeModel\">\n" +
            "    <tree idref=\"startingTree\"/>\n" +
            "    <rootHeight><parameter id=\"treeModel.rootHeight\"/></rootHeight>\n" +
            "    <nodeHeights internalNodes=\"true\"><parameter id=\"treeModel.internalNodeHeights\"/></nodeHeights>\n" +
            "    <nodeHeights internalNodes=\"true\" rootNode=\"true\">" +
            "<parameter id=\"treeModel.allInternalNodeHeights\"/></nodeHeights>\n" +
            "  </treeModel>\n" +
            "  <tmrcaStatistic id=\"tmrca(clade)\"><treeModel idref=\"treeModel\"/>" +
            "<mrca><taxa idref=\"clade\"/></mrca></tmrcaStatistic>\n";

    public TMRCAPriorGradientRefusalTest(String name) {
        super(name);
    }

    public void testGradientOfTMRCAPrior() {
        String xml = "<beast>\n" + TREE +
                "  <logNormalPrior id=\"tmrcaPrior\" mu=\"0.0\" sigma=\"1.0\" offset=\"0.0\">" +
                "<statistic idref=\"tmrca(clade)\"/></logNormalPrior>\n" +
                "  <gradient id=\"gradient\"><logNormalPrior idref=\"tmrcaPrior\"/>" +
                "<parameter idref=\"treeModel.allInternalNodeHeights\"/></gradient>\n" +
                "</beast>\n";
        assertRefused(xml, "tMRCA calibration priors");
    }

    public void testSpeciationGradientWithCalibration() {
        String xml = "<beast>\n" + TREE +
                "  <yuleModel id=\"yule\" units=\"substitutions\">" +
                "<birthRate><parameter id=\"yule.birthRate\" value=\"1.0\"/></birthRate></yuleModel>\n" +
                "  <speciationLikelihood id=\"speciation\">\n" +
                "    <model><yuleModel idref=\"yule\"/></model>\n" +
                "    <speciesTree><treeModel idref=\"treeModel\"/></speciesTree>\n" +
                "    <calibration><point><uniformDistributionModel><lower>0.5</lower><upper>1.5</upper>" +
                "</uniformDistributionModel><taxa idref=\"clade\"/></point></calibration>\n" +
                "  </speciationLikelihood>\n" +
                "  <speciationLikelihoodGradient id=\"gradient\" wrtParameter=\"nodeHeight\">" +
                "<speciationLikelihood idref=\"speciation\"/><treeModel idref=\"treeModel\"/>" +
                "</speciationLikelihoodGradient>\n" +
                "</beast>\n";
        assertRefused(xml, "tMRCA calibration priors");
    }

    private static void assertRefused(String xml, String expectedMessage) {
        BeastParser parser = new BeastParser(new String[]{}, null, false, false, false, new BeastVersion());
        try {
            parser.parse(new StringReader(xml), null);
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
