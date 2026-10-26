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
import dr.evomodel.substmodel.RelativeRatesProvider;
import dr.evomodel.substmodel.nucleotide.HKY;
import dr.evomodel.substmodel.nucleotide.TN93;
import dr.evomodelxml.substmodel.FrequencyModelParser;
import dr.evomodelxml.substmodel.HKYParser;
import dr.evomodelxml.substmodel.HierarchicalSubstitutionModelParser;
import dr.evomodelxml.substmodel.TN93Parser;
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
        parser.addXMLObjectParser(new TN93Parser());
        parser.addXMLObjectParser(new dr.evomodelxml.substmodel.BinarySubstitutionModelParser());
        parser.addXMLObjectParser(new dr.evomodelxml.substmodel.GeneralSubstitutionModelParser());
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
        assertTrue(message, message.contains("characterTree is required with ratesParameter on data type nucleotide"));
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

    private static final String SHARED_FREQUENCIES = "    <frequencies><frequencyModel idref=\"hky.frequencies.model\"/></frequencies>\n";
    private static final String HKY_SHARED =
            "  <HKYModel id=\"hky\">\n" +
            "    <frequencies>\n" +
            "      <frequencyModel id=\"hky.frequencies.model\" dataType=\"nucleotide\">\n" +
            "        <frequencies><parameter id=\"hky.frequencies\" value=\"0.1 0.2 0.3 0.4\"/></frequencies>\n" +
            "      </frequencyModel>\n" +
            "    </frequencies>\n" +
            "    <kappa><parameter id=\"kappa\" value=\"2.5\" lower=\"0.0\"/></kappa>\n" +
            "  </HKYModel>\n";

    private static String hrate(String attributes, String body) {
        return "  <hierarchicalSubstitutionModel id=\"hrate\"" + attributes + ">\n" + body + "  </hierarchicalSubstitutionModel>\n";
    }

    public void testKappaForm() throws Exception {
        for (String attributes : new String[]{"", " characterTreeNewick=\"((T:1,C:1):1,(G:1,A:1):1);\""}) {
            Map<String, XMLObject> store = parse(HKY_SHARED + hrate(attributes,
                    SHARED_FREQUENCIES + "    <kappa><parameter idref=\"kappa\"/></kappa>\n"));
            HierarchicalSubstitutionModel model = (HierarchicalSubstitutionModel) store.get("hrate").getNativeObject();
            HKY hky = (HKY) store.get("hky").getNativeObject();
            assertTrue(model.getRateProvider() instanceof RelativeRatesProvider);
            RelativeRatesProvider provider = (RelativeRatesProvider) model.getRateProvider();
            assertEquals(3, provider.getCladeNames().length);
            assertEquals(0.5825242718446603, provider.getLambda("AG"), 1e-15);
            assertEquals(0.0, maxPDiff(model, hky), 1e-13);
        }
    }

    public void testKappaWithCharacterTree() throws Exception {
        Map<String, XMLObject> store = parse(TAXA_AND_TREE + HKY_SHARED + hrate("",
                SHARED_FREQUENCIES + CHARACTER_TREE + "    <kappa><parameter idref=\"kappa\"/></kappa>\n"));
        HierarchicalSubstitutionModel model = (HierarchicalSubstitutionModel) store.get("hrate").getNativeObject();
        assertEquals(0.0, maxPDiff(model, (HKY) store.get("hky").getNativeObject()), 1e-13);
    }

    public void testTN93Form() throws Exception {
        Map<String, XMLObject> store = parse(
                "  <TN93Model id=\"tn93\">\n" +
                "    <frequencies><frequencyModel id=\"f\" dataType=\"nucleotide\"><frequencies><parameter value=\"0.1 0.2 0.3 0.4\"/></frequencies></frequencyModel></frequencies>\n" +
                "    <kappa1><parameter id=\"k1\" value=\"1.7\"/></kappa1>\n" +
                "    <kappa2><parameter id=\"k2\" value=\"3.1\"/></kappa2>\n" +
                "  </TN93Model>\n" +
                hrate("", "    <frequencies><frequencyModel idref=\"f\"/></frequencies>\n" +
                        "    <kappa1><parameter idref=\"k1\"/></kappa1><kappa2><parameter idref=\"k2\"/></kappa2>\n"));
        HierarchicalSubstitutionModel model = (HierarchicalSubstitutionModel) store.get("hrate").getNativeObject();
        TN93 tn93 = (TN93) store.get("tn93").getNativeObject();
        double[] ph = new double[16], pd = new double[16];
        for (double t : new double[]{0.05, 0.37, 1.9}) {
            model.getTransitionProbabilities(t, ph);
            tn93.getTransitionProbabilities(t, pd);
            for (int i = 0; i < 16; i++)
                assertEquals(pd[i], ph[i], 1e-13);
        }
    }

    public void testSelectorErrors() throws Exception {
        String message = parseFailure(TAXA_AND_TREE + hrate("", FREQUENCIES + CHARACTER_TREE + RATES
                + "    <kappa><parameter value=\"2.0\"/></kappa>\n"));
        assertTrue(message, message.contains("give at most one of ratesParameter, kappa"));

        message = parseFailure(hrate("", FREQUENCIES + "    <kappa1><parameter value=\"2.0\"/></kappa1>\n"));
        assertTrue(message, message.contains("kappa1 and kappa2 must be given together"));

        message = parseFailure(hrate("",
                "    <frequencies><frequencyModel dataType=\"binary\"><frequencies><parameter value=\"0.3 0.7\"/></frequencies></frequencyModel></frequencies>\n"
                        + "    <kappa><parameter value=\"2.0\"/></kappa>\n"));
        assertTrue(message, message.contains("a nucleotide frequency model"));

        message = parseFailure(hrate(" characterTreeNewick=\"((A:1,C:1):1,(G:1,T:1):1);\"",
                FREQUENCIES + "    <kappa><parameter value=\"1.0\"/></kappa>\n"));
        assertTrue(message, message.contains("the character tree must have the clades {A,G} and {C,T}"));

        message = parseFailure(TAXA_AND_TREE + hrate(" characterTreeNewick=\"((A:1,G:1):1,(C:1,T:1):1);\"",
                FREQUENCIES + CHARACTER_TREE + "    <kappa><parameter value=\"2.0\"/></kappa>\n"));
        assertTrue(message, message.contains("characterTree or characterTreeNewick, not both"));

        message = parseFailure(hrate("", FREQUENCIES));
        assertTrue(message, message.contains("characterTree is required with no rate selector"));
    }

    public void testRatesFromHKY() throws Exception {
        for (String frequencies : new String[]{"", SHARED_FREQUENCIES}) {
            Map<String, XMLObject> store = parse(HKY_SHARED + hrate("",
                    frequencies + "    <ratesFrom><HKYModel idref=\"hky\"/></ratesFrom>\n"));
            HierarchicalSubstitutionModel model = (HierarchicalSubstitutionModel) store.get("hrate").getNativeObject();
            HKY hky = (HKY) store.get("hky").getNativeObject();
            assertSame(hky.getFrequencyModel(), model.getFrequencyModel());
            RelativeRatesProvider provider = (RelativeRatesProvider) model.getRateProvider();
            assertEquals(0.5825242718446603, provider.getLambda("AG"), 1e-15);
            assertEquals(0.0, maxPDiff(model, hky), 1e-13);
            ((dr.inference.model.Parameter) store.get("kappa").getNativeObject()).setParameterValue(0, 7.0);
            assertEquals(0.0, maxPDiff(model, hky), 1e-13);
        }
    }

    public void testRatesFromTN93() throws Exception {
        Map<String, XMLObject> store = parse(
                "  <TN93Model id=\"tn93\">\n" +
                "    <frequencies><frequencyModel id=\"f\" dataType=\"nucleotide\"><frequencies><parameter value=\"0.1 0.2 0.3 0.4\"/></frequencies></frequencyModel></frequencies>\n" +
                "    <kappa1><parameter id=\"k1\" value=\"1.7\"/></kappa1>\n" +
                "    <kappa2><parameter id=\"k2\" value=\"3.1\"/></kappa2>\n" +
                "  </TN93Model>\n" +
                hrate(" characterTreeNewick=\"((T:1,C:1):1,(G:1,A:1):1);\"", "    <ratesFrom><TN93Model idref=\"tn93\"/></ratesFrom>\n"));
        HierarchicalSubstitutionModel model = (HierarchicalSubstitutionModel) store.get("hrate").getNativeObject();
        TN93 tn93 = (TN93) store.get("tn93").getNativeObject();
        double[] ph = new double[16], pd = new double[16];
        for (double t : new double[]{0.05, 0.37, 1.9}) {
            model.getTransitionProbabilities(t, ph);
            tn93.getTransitionProbabilities(t, pd);
            for (int i = 0; i < 16; i++)
                assertEquals(pd[i], ph[i], 1e-13);
        }
    }

    public void testRatesFromErrors() throws Exception {
        String message = parseFailure(HKY_SHARED + hrate("", FREQUENCIES + "    <ratesFrom><HKYModel idref=\"hky\"/></ratesFrom>\n"));
        assertTrue(message, message.contains("the frequency model of ratesFrom's model must be"));

        message = parseFailure(HKY_SHARED + hrate("", "    <ratesFrom><HKYModel idref=\"hky\"/></ratesFrom>\n"
                + "    <kappa><parameter value=\"2.0\"/></kappa>\n"));
        assertTrue(message, message.contains("give at most one of") && message.contains("ratesFrom"));

        message = parseFailure(HKY_SHARED + hrate(" characterTreeNewick=\"((A:1,C:1):1,(G:1,T:1):1);\"",
                "    <ratesFrom><HKYModel idref=\"hky\"/></ratesFrom>\n"));
        assertTrue(message, message.contains("the character tree must have the clades {A,G} and {C,T}"));

        message = parseFailure(hrate("", ""));
        assertTrue(message, message.contains("frequencies is required unless ratesFrom gives the model"));
    }

    private static double maxPDiff(dr.evomodel.substmodel.SubstitutionModel a, dr.evomodel.substmodel.SubstitutionModel b) {
        final int n = a.getDataType().getStateCount();
        double[] pa = new double[n * n], pb = new double[n * n];
        double d = 0.0;
        for (double t : new double[]{0.0, 0.05, 0.37, 1.9}) {
            a.getTransitionProbabilities(t, pa);
            b.getTransitionProbabilities(t, pb);
            for (int i = 0; i < n * n; i++)
                d = Math.max(d, Math.abs(pa[i] - pb[i]));
        }
        return d;
    }

    public void testEqualRatesNucleotides() throws Exception {
        Map<String, XMLObject> store = parse(
                "  <HKYModel id=\"f81\">\n" +
                "    <frequencies><frequencyModel id=\"f\" dataType=\"nucleotide\"><frequencies><parameter value=\"0.1 0.2 0.3 0.4\"/></frequencies></frequencyModel></frequencies>\n" +
                "    <kappa><parameter value=\"1.0\"/></kappa>\n" +
                "  </HKYModel>\n" +
                hrate("", "    <frequencies><frequencyModel idref=\"f\"/></frequencies>\n    <equalRates/>\n"));
        HierarchicalSubstitutionModel model = (HierarchicalSubstitutionModel) store.get("hrate").getNativeObject();
        assertTrue(model.getRateProvider() instanceof dr.evomodel.substmodel.EqualRatesProvider);
        assertEquals(0.0, maxPDiff(model, (HKY) store.get("f81").getNativeObject()), 1e-13);
        assertEquals(1.0 / (1.0 - 0.3), model.getBeagleModelLambdas()[3], 1e-14);
    }

    private static final String BINARY =
            "  <binarySubstitutionModel id=\"bin\">\n" +
            "    <frequencies><frequencyModel id=\"bf\" dataType=\"binary\"><frequencies><parameter value=\"0.3 0.7\"/></frequencies></frequencyModel></frequencies>\n" +
            "  </binarySubstitutionModel>\n";

    public void testEqualRatesBinary() throws Exception {
        Map<String, XMLObject> store = parse(BINARY
                + hrate("", "    <frequencies><frequencyModel idref=\"bf\"/></frequencies>\n    <equalRates/>\n")
                .replace("id=\"hrate\"", "id=\"hrate.equal\"")
                + hrate("", "    <ratesFrom><binarySubstitutionModel idref=\"bin\"/></ratesFrom>\n")
                .replace("id=\"hrate\"", "id=\"hrate.from\""));
        dr.evomodel.substmodel.SubstitutionModel bin = (dr.evomodel.substmodel.SubstitutionModel) store.get("bin").getNativeObject();
        for (String id : new String[]{"hrate.equal", "hrate.from"}) {
            HierarchicalSubstitutionModel model = (HierarchicalSubstitutionModel) store.get(id).getNativeObject();
            assertEquals(id, 2, model.getDataType().getStateCount());
            assertEquals(id, 0.0, maxPDiff(model, bin), 1e-13);
            assertEquals(id, 1.0 / (2 * 0.3 * 0.7), model.getBeagleModelLambdas()[1], 1e-13);
        }
    }

    public void testEqualRatesNotNormalised() throws Exception {
        Map<String, XMLObject> store = parse(hrate("", FREQUENCIES + "    <equalRates normalised=\"false\"/>\n"));
        HierarchicalSubstitutionModel model = (HierarchicalSubstitutionModel) store.get("hrate").getNativeObject();
        assertEquals(1.0, model.getBeagleModelLambdas()[3], 0.0);
    }

    public void testEqualRatesErrors() throws Exception {
        String message = parseFailure(hrate("", FREQUENCIES + "    <equalRates/>\n    <kappa><parameter value=\"2.0\"/></kappa>\n"));
        assertTrue(message, message.contains("give at most one of") && message.contains("equalRates"));

        message = parseFailure(
                "  <generalSubstitutionModel id=\"gsm\">\n" +
                "    <frequencies><frequencyModel id=\"gf\" dataType=\"nucleotide\"><frequencies><parameter value=\"0.1 0.2 0.3 0.4\"/></frequencies></frequencyModel></frequencies>\n" +
                "    <rates><parameter value=\"1 1 1 1 1 2\"/></rates>\n" +
                "  </generalSubstitutionModel>\n" +
                hrate("", "    <ratesFrom><generalSubstitutionModel idref=\"gsm\"/></ratesFrom>\n"));
        assertTrue(message, message.contains("does not have equal rates"));
    }
}
