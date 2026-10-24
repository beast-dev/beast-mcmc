/*
 * HierarchicalSubstitutionModelParserTest.java
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

package test.dr.evomodel.substmodel;

import dr.evomodel.substmodel.HierarchicalSubstitutionModel;
import dr.evomodel.substmodel.nucleotide.HKY;
import dr.evomodelxml.substmodel.FrequencyModelParser;
import dr.evomodelxml.substmodel.HKYParser;
import dr.evomodelxml.substmodel.HierarchicalSubstitutionModelParser;
import dr.evomodelxml.tree.TreeModelParser;
import dr.evoxml.NewickParser;
import dr.evoxml.TaxaParser;
import dr.evoxml.TaxonParser;
import dr.inference.model.ParameterParser;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLParser;
import junit.framework.TestCase;

import java.io.StringReader;
import java.util.Map;

/**
 * Parser round trip and error paths of hierarchicalSubstitutionModel.
 *
 * @author Xiang Ji
 */
public class HierarchicalSubstitutionModelParserTest extends TestCase {

    private static final String TAXA_AND_TREE =
            "  <taxa id=\"characterTaxa\">\n" +
            "    <taxon id=\"A\"/><taxon id=\"C\"/><taxon id=\"G\"/><taxon id=\"T\"/>\n" +
            "  </taxa>\n" +
            "  <newick id=\"characterNewick\" usingHeights=\"true\" usingDates=\"false\">\n" +
            "    ((A:1,G:1):1,(C:1,T:1):1);\n" +
            "  </newick>\n" +
            "  <treeModel id=\"characterTree\" fixHeights=\"true\">\n" +
            "    <tree idref=\"characterNewick\"/>\n" +
            "    <rootHeight><parameter id=\"characterTree.rootHeight\"/></rootHeight>\n" +
            "    <nodeHeights internalNodes=\"true\"><parameter id=\"characterTree.internalNodeHeights\"/></nodeHeights>\n" +
            "  </treeModel>\n";

    private static final String FREQUENCIES =
            "    <frequencies>\n" +
            "      <frequencyModel dataType=\"nucleotide\">\n" +
            "        <frequencies><parameter id=\"hrate.frequencies\" value=\"0.1 0.2 0.3 0.4\"/></frequencies>\n" +
            "      </frequencyModel>\n" +
            "    </frequencies>\n";

    private static final String HKY_MODEL =
            "  <HKYModel id=\"hky\">\n" +
            "    <frequencies>\n" +
            "      <frequencyModel dataType=\"nucleotide\">\n" +
            "        <frequencies><parameter id=\"hky.frequencies\" value=\"0.1 0.2 0.3 0.4\"/></frequencies>\n" +
            "      </frequencyModel>\n" +
            "    </frequencies>\n" +
            "    <kappa><parameter id=\"kappa\" value=\"2.5\" lower=\"0.0\"/></kappa>\n" +
            "  </HKYModel>\n";

    private static final String CHARACTER_TREE = "    <characterTree><treeModel idref=\"characterTree\"/></characterTree>\n";
    private static final String RATES =
            "    <ratesParameter><parameter id=\"hrate.rates\" value=\"0.5825242718446603 0.8737864077669906 0.9708737864077672\"/></ratesParameter>\n";
    private static final String SCALED_RATES =
            "    <ratesParameter normalised=\"true\"><parameter id=\"hrate.rates\" value=\"2.155339805825243 3.2330097087378652 3.592233009708739\"/></ratesParameter>\n";

    private static String model(String tree, String rates) {
        return "  <hierarchicalSubstitutionModel id=\"hrate\">\n" + FREQUENCIES + tree + rates + "  </hierarchicalSubstitutionModel>\n";
    }

    private static Map<String, XMLObject> parse(String body) throws Exception {
        XMLParser parser = new XMLParser(true, true, false, null);
        parser.addXMLObjectParser(new TaxaParser());
        parser.addXMLObjectParser(new TaxonParser());
        parser.addXMLObjectParser(new NewickParser());
        parser.addXMLObjectParser(new TreeModelParser());
        parser.addXMLObjectParser(new ParameterParser());
        parser.addXMLObjectParser(new FrequencyModelParser());
        parser.addXMLObjectParser(new HKYParser());
        parser.addXMLObjectParser(new HierarchicalSubstitutionModelParser());
        return parser.parse(new StringReader("<?xml version=\"1.0\" standalone=\"yes\"?>\n<beast>\n" + body + "</beast>\n"), false);
    }

    private static String parseFailure(String body) throws Exception {
        try {
            parse(body);
        } catch (XMLParseException e) {
            return e.getMessage();
        }
        fail("the parse should have failed with an XMLParseException");
        return null;
    }

    private static double maxPDiff(HierarchicalSubstitutionModel hrate, HKY hky) {
        double[] ph = new double[16], pd = new double[16];
        double d = 0.0;
        for (double t : new double[]{0.0, 0.05, 0.37, 1.9}) {
            hrate.getTransitionProbabilities(t, ph);
            hky.getTransitionProbabilities(t, pd);
            for (int i = 0; i < 16; i++)
                d = Math.max(d, Math.abs(ph[i] - pd[i]));
        }
        return d;
    }

    public void testRoundTrip() throws Exception {
        Map<String, XMLObject> store = parse(TAXA_AND_TREE + model(CHARACTER_TREE, RATES) + HKY_MODEL);
        HierarchicalSubstitutionModel hrate = (HierarchicalSubstitutionModel) store.get("hrate").getNativeObject();
        HKY hky = (HKY) store.get("hky").getNativeObject();
        assertTrue(hrate.getRateProvider() instanceof HierarchicalSubstitutionModel.HierarchicalRateProvider.Default);
        assertEquals(0.0, maxPDiff(hrate, hky), 1e-13);
        double[] lambdas = hrate.getBeagleModelLambdas();
        assertEquals(65539.0, lambdas[0]);
        assertEquals(0.5825242718446603, lambdas[1], 1e-15);
        assertEquals(0.8737864077669906, lambdas[2], 1e-15);
        assertEquals(0.9708737864077672, lambdas[3], 1e-15);
    }

    public void testNormalisedRates() throws Exception {
        Map<String, XMLObject> store = parse(TAXA_AND_TREE + model(CHARACTER_TREE, SCALED_RATES) + HKY_MODEL);
        HierarchicalSubstitutionModel hrate = (HierarchicalSubstitutionModel) store.get("hrate").getNativeObject();
        HKY hky = (HKY) store.get("hky").getNativeObject();
        assertEquals(0.0, maxPDiff(hrate, hky), 1e-13);
    }

    public void testMissingCharacterTree() throws Exception {
        String message = parseFailure(TAXA_AND_TREE + model("", RATES));
        assertTrue(message, message.contains("hierarchicalSubstitutionModel") && message.contains("characterTree"));
    }

    public void testMissingFrequencies() throws Exception {
        String message = parseFailure(TAXA_AND_TREE +
                "  <hierarchicalSubstitutionModel id=\"hrate\">\n" + CHARACTER_TREE + RATES + "  </hierarchicalSubstitutionModel>\n");
        assertTrue(message, message.contains("hierarchicalSubstitutionModel") && message.contains("frequencies"));
    }

    public void testWrongRatesDimension() throws Exception {
        String message = parseFailure(TAXA_AND_TREE + model(CHARACTER_TREE,
                "    <ratesParameter><parameter id=\"hrate.rates\" value=\"0.5 0.9\"/></ratesParameter>\n"));
        assertTrue(message, message.contains("ratesParameter has dimension 2 but the character tree has 3 internal nodes"));
    }

    public void testTipNotAState() throws Exception {
        String message = parseFailure(TAXA_AND_TREE.replace("<taxon id=\"T\"/>", "<taxon id=\"X\"/>").replace("T:1", "X:1")
                + model(CHARACTER_TREE, RATES));
        assertTrue(message, message.contains("Character tree tip X is not a distinct state of nucleotide"));
    }

    public void testWrongDataType() throws Exception {
        String message = parseFailure(TAXA_AND_TREE +
                "  <hierarchicalSubstitutionModel id=\"hrate\">\n" +
                "    <frequencies><frequencyModel dataType=\"binary\"><frequencies><parameter value=\"0.3 0.7\"/></frequencies></frequencyModel></frequencies>\n" +
                CHARACTER_TREE + RATES + "  </hierarchicalSubstitutionModel>\n");
        assertTrue(message, message.contains("Character tree has 4 tips but the data type has 2 states"));
    }
}
