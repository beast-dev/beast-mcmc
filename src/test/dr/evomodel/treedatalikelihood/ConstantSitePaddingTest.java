/*
 * ConstantSitePaddingTest.java
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
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.Likelihood;
import dr.xml.XMLObject;
import junit.framework.TestCase;

import java.io.StringReader;
import java.util.Map;

/**
 * Constant-site padding (plan 8 M4a, IQ-TREE -fconst): the variable columns plus counts of constant columns, merged
 * with <constantPatterns> and <mergePatterns>, give the same likelihood and analytic gradient as the full alignment,
 * since the counts are pattern weights and the constant columns' likelihoods do not depend on which columns they
 * were. And treeDataLikelihood refuses ascertained patterns, whose correction it would ignore (M4-guard).
 */
public class ConstantSitePaddingTest extends TestCase {

    private static final String[] TAXA = {"a", "b", "c", "d", "e"};
    // the variable columns, then constant columns: 4 A, 3 C, 2 G and 5 T
    private static final String[] VARIABLE = {"ACGTTGCAAG", "ACGTTGCATG", "ACCTTGGATG", "TCCTAGGATC", "TCCAAGGCTC"};
    private static final String CONSTANT = "AAAACCCGGTTTTT";
    private static final String COUNTS = "4 3 2 5";

    public ConstantSitePaddingTest(String name) {
        super(name);
    }

    public void testPaddingEqualsTheFullAlignment() throws Exception {
        Map<String, XMLObject> full = parse(alignment(true), "<patterns id=\"patterns\"><alignment idref=\"alignment\"/></patterns>");
        Map<String, XMLObject> padded = parse(alignment(false),
                "<mergePatterns id=\"patterns\">" +
                        "<patterns id=\"variable\"><alignment idref=\"alignment\"/></patterns>" +
                        "<constantPatterns><patterns idref=\"variable\"/>" +
                        "<counts><parameter value=\"" + COUNTS + "\"/></counts></constantPatterns>" +
                        "</mergePatterns>");

        double fullLogL = logL(full);
        double paddedLogL = logL(padded);
        assertEquals(fullLogL, paddedLogL, 1E-8 * Math.abs(fullLogL));
        // without the padding the likelihood is that of the variable columns alone
        Map<String, XMLObject> unpadded = parse(alignment(false), "<patterns id=\"patterns\"><alignment idref=\"alignment\"/></patterns>");
        assertTrue(Math.abs(logL(unpadded) - fullLogL) > 1.0);

        double[] fullGradient = gradient(full);
        double[] paddedGradient = gradient(padded);
        assertEquals(fullGradient[0], paddedGradient[0], 1E-8 * Math.max(1.0, Math.abs(fullGradient[0])));
    }

    public void testAscertainedPatternsAreRefused() {
        try {
            parse(alignment(false), "<ascertainedPatterns id=\"patterns\"><alignment idref=\"alignment\"/>" +
                    "<excludePatterns from=\"1\" to=\"4\"/></ascertainedPatterns>");
            fail("expected treeDataLikelihood to refuse ascertained patterns");
        } catch (Exception e) {
            String message = e.toString();
            for (Throwable cause = e.getCause(); cause != null; cause = cause.getCause()) {
                message += " / " + cause;
            }
            assertTrue(message, message.contains("ascertainment correction") && message.contains("constantPatterns"));
        }
    }

    private static String alignment(boolean withConstantColumns) {
        StringBuilder sb = new StringBuilder("<alignment id=\"alignment\" dataType=\"nucleotide\">\n");
        for (int i = 0; i < TAXA.length; i++) {
            sb.append("  <sequence><taxon idref=\"").append(TAXA[i]).append("\"/>").append(VARIABLE[i])
                    .append(withConstantColumns ? CONSTANT : "").append("</sequence>\n");
        }
        return sb.append("</alignment>\n").toString();
    }

    private static Map<String, XMLObject> parse(String alignment, String patterns) throws Exception {
        StringBuilder taxa = new StringBuilder("<taxa id=\"taxa\">");
        for (String taxon : TAXA) {
            taxa.append("<taxon id=\"").append(taxon).append("\"/>");
        }
        String xml = "<beast>\n" + taxa + "</taxa>\n" + alignment + patterns + "\n" +
                "<newick id=\"tree\">((a:0.1,b:0.2):0.15,(c:0.3,(d:0.1,e:0.25):0.05):0.1);</newick>\n" +
                "<treeModel id=\"treeModel\"><tree idref=\"tree\"/>" +
                "<rootHeight><parameter id=\"rootHeight\"/></rootHeight>" +
                "<nodeHeights internalNodes=\"true\"><parameter id=\"internalHeights\"/></nodeHeights>" +
                "<nodeHeights internalNodes=\"true\" rootNode=\"true\"><parameter id=\"allHeights\"/></nodeHeights>" +
                "</treeModel>\n" +
                "<strictClockBranchRates id=\"branchRates\"><rate><parameter id=\"clock.rate\" value=\"1.0\"/></rate>" +
                "</strictClockBranchRates>\n" +
                "<HKYModel id=\"hky\"><frequencies><frequencyModel dataType=\"nucleotide\">" +
                "<frequencies><parameter id=\"frequencies\" value=\"0.3 0.2 0.2 0.3\"/></frequencies></frequencyModel>" +
                "</frequencies><kappa><parameter id=\"kappa\" value=\"3.0\" lower=\"0.0\"/></kappa></HKYModel>\n" +
                "<siteModel id=\"siteModel\"><substitutionModel><HKYModel idref=\"hky\"/></substitutionModel>" +
                "<gammaShape gammaCategories=\"4\"><parameter id=\"alpha\" value=\"0.7\"/></gammaShape></siteModel>\n" +
                "<treeDataLikelihood id=\"treeLikelihood\" usePreOrder=\"true\" branchInfinitesimalDerivative=\"true\">" +
                "<partition><patterns idref=\"patterns\"/><siteModel idref=\"siteModel\"/></partition>" +
                "<treeModel idref=\"treeModel\"/><strictClockBranchRates idref=\"branchRates\"/></treeDataLikelihood>\n" +
                "<branchSubstitutionParameterGradient id=\"kappaGradient\" homogeneous=\"true\">" +
                "<treeDataLikelihood idref=\"treeLikelihood\"/><parameter idref=\"kappa\"/>" +
                "</branchSubstitutionParameterGradient>\n" +
                "</beast>\n";
        xml = xml.replace("<patterns idref=\"patterns\"/>", patterns.startsWith("<mergePatterns") ?
                "<mergePatterns idref=\"patterns\"/>" : patterns.startsWith("<ascertainedPatterns") ?
                "<ascertainedPatterns idref=\"patterns\"/>" : "<patterns idref=\"patterns\"/>");
        BeastParser parser = new BeastParser(new String[]{}, null, false, false, false, new BeastVersion());
        return parser.parse(new StringReader(xml), false);
    }

    private static double logL(Map<String, XMLObject> objects) {
        Likelihood likelihood = (Likelihood) objects.get("treeLikelihood").getNativeObject();
        likelihood.makeDirty();
        return likelihood.getLogLikelihood();
    }

    private static double[] gradient(Map<String, XMLObject> objects) {
        return ((GradientWrtParameterProvider) objects.get("kappaGradient").getNativeObject()).getGradientLogDensity();
    }
}
