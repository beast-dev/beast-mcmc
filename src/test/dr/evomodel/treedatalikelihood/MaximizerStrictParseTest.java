/*
 * MaximizerStrictParseTest.java
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

package test.dr.evomodel.treedatalikelihood;

import dr.app.beast.BeastParser;
import dr.app.beast.BeastVersion;
import dr.evomodel.treedatalikelihood.discrete.MaximizerWrtParameter;
import dr.xml.XMLObject;
import junit.framework.TestCase;

import java.io.StringReader;
import java.util.Map;

/**
 * Under -strict, <maximizeWrtParameter> accepts any multivariable transform, including <nodeHeightTransform>, and
 * <treeDataLikelihood> accepts the pre-order attributes that gradients need; maximizeOnParse="false" parses the
 * maximizer without running it.
 */
public class MaximizerStrictParseTest extends TestCase {

    private static final String XML = "<?xml version=\"1.0\" standalone=\"yes\"?>\n" +
            "<beast>\n" +
            "  <taxa id=\"taxa\"><taxon id=\"A\"/><taxon id=\"B\"/><taxon id=\"C\"/><taxon id=\"D\"/></taxa>\n" +
            "  <alignment id=\"alignment\" dataType=\"nucleotide\">\n" +
            "    <sequence><taxon idref=\"A\"/>ACGTACGTTACGATCGATCG</sequence>\n" +
            "    <sequence><taxon idref=\"B\"/>ACGTACGTAACGATCGTTCG</sequence>\n" +
            "    <sequence><taxon idref=\"C\"/>ACGAACGTTACGTTCGATCG</sequence>\n" +
            "    <sequence><taxon idref=\"D\"/>TCGTACGATACGATCGATCC</sequence>\n" +
            "  </alignment>\n" +
            "  <patterns id=\"patterns\" from=\"1\" strip=\"false\"><alignment idref=\"alignment\"/></patterns>\n" +
            "  <newick id=\"startingTree\" usingDates=\"false\" usingHeights=\"true\">" +
            "((A:0.1,B:0.1):0.05,(C:0.07,D:0.07):0.08);</newick>\n" +
            "  <treeModel id=\"treeModel\">\n" +
            "    <tree idref=\"startingTree\"/>\n" +
            "    <rootHeight><parameter id=\"treeModel.rootHeight\"/></rootHeight>\n" +
            "    <nodeHeights internalNodes=\"true\"><parameter id=\"treeModel.internalNodeHeights\"/></nodeHeights>\n" +
            "    <nodeHeights internalNodes=\"true\" rootNode=\"true\">" +
            "<parameter id=\"treeModel.allInternalNodeHeights\"/></nodeHeights>\n" +
            "  </treeModel>\n" +
            "  <strictClockBranchRates id=\"branchRates\"><rate><parameter id=\"clock.rate\" value=\"1.0\"/></rate>" +
            "</strictClockBranchRates>\n" +
            "  <HKYModel id=\"hky\">\n" +
            "    <frequencies><frequencyModel dataType=\"nucleotide\">\n" +
            "      <frequencies><parameter id=\"frequencies\" value=\"0.25 0.25 0.25 0.25\"/></frequencies>\n" +
            "    </frequencyModel></frequencies>\n" +
            "    <kappa><parameter id=\"kappa\" value=\"2.0\" lower=\"0.0\"/></kappa>\n" +
            "  </HKYModel>\n" +
            "  <siteModel id=\"siteModel\"><substitutionModel><HKYModel idref=\"hky\"/></substitutionModel></siteModel>\n" +
            "  <treeDataLikelihood id=\"treeLikelihood\" useAmbiguities=\"false\" usePreOrder=\"true\"" +
            " branchRateDerivative=\"true\" branchInfinitesimalDerivative=\"true\">\n" +
            "    <partition><patterns idref=\"patterns\"/><siteModel idref=\"siteModel\"/></partition>\n" +
            "    <treeModel idref=\"treeModel\"/>\n" +
            "    <strictClockBranchRates idref=\"branchRates\"/>\n" +
            "  </treeDataLikelihood>\n" +
            "  <nodeHeightTransform id=\"nodeHeightTransform\" withRoot=\"true\">\n" +
            "    <tree idref=\"treeModel\"/>\n" +
            "    <strictClockBranchRates idref=\"branchRates\"/>\n" +
            "    <nodeHeights><nodeHeightProxyParameter id=\"nodeHeights\" rootNode=\"true\">" +
            "<tree idref=\"treeModel\"/></nodeHeightProxyParameter></nodeHeights>\n" +
            "    <ratios realLine=\"true\"><parameter id=\"nodeHeight.ratios\"/></ratios>\n" +
            "  </nodeHeightTransform>\n" +
            "  <maximizeWrtParameter id=\"maximizer\" nIterations=\"20\" maximizeOnParse=\"false\">\n" +
            "    <densityWrtParameter>\n" +
            "      <treeDataLikelihood idref=\"treeLikelihood\"/>\n" +
            "      <nodeHeightProxyParameter idref=\"nodeHeights\"/>\n" +
            "    </densityWrtParameter>\n" +
            "    <nodeHeightTransform idref=\"nodeHeightTransform\"/>\n" +
            "  </maximizeWrtParameter>\n" +
            "</beast>\n";

    public MaximizerStrictParseTest(String name) {
        super(name);
    }

    public void testStrictParseWithoutMaximizing() throws Exception {
        BeastParser parser = new BeastParser(new String[]{}, null, false, false, true, new BeastVersion());
        Map<String, XMLObject> objects = parser.parse(new StringReader(XML), false);
        MaximizerWrtParameter maximizer = (MaximizerWrtParameter) objects.get("maximizer").getNativeObject();

        assertFalse(maximizer.wasExecuted());
        assertEquals("Not yet executed.", maximizer.getReport());

        final double before = maximizer.getLikelihood().getLogLikelihood();
        maximizer.maximize();

        assertTrue(maximizer.wasExecuted());
        assertNotNull(maximizer.getReason());
        assertTrue(maximizer.getLikelihood().getLogLikelihood() >= before);
    }
}
