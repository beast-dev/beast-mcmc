/*
 * GeneralOUTraitSimulator.java
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
 *
 */

package dr.evomodel.continuous;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxon;
import dr.inferencexml.loggers.LoggerParser;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.AttributeRule;
import dr.xml.StringAttributeRule;
import dr.xml.XMLObject;
import dr.xml.XMLObjectParser;
import dr.xml.XMLParseException;
import dr.xml.XMLParser;
import dr.xml.XMLSyntaxRule;

import java.io.PrintWriter;

/**
 * XML-facing exact OU simulator that writes simulated tip traits as a reusable
 * BEAST taxa block.
 *
 * <p>The actual OU draw is delegated to {@link GeneralOUTreeSimulator}, so this
 * class uses the same canonical transition kernel as the OU likelihood.</p>
 */
public final class GeneralOUTraitSimulator {

    public static final String GENERAL_OU_TRAIT_SIMULATOR = "generalOuTraitSimulator";
    public static final String TAXA_ID = "taxaId";
    public static final String FILE_NAME = "fileName";

    private GeneralOUTraitSimulator() {
        // Utility class.
    }

    private static void writeTaxaBlock(final Tree tree,
                                       final String traitName,
                                       final String taxaId,
                                       final PrintWriter writer) throws XMLParseException {
        writer.println("<taxa id=\"" + escapeXml(taxaId) + "\">");
        for (int i = 0; i < tree.getExternalNodeCount(); ++i) {
            final NodeRef node = tree.getExternalNode(i);
            final Taxon taxon = tree.getNodeTaxon(node);
            if (taxon == null) {
                throw new XMLParseException("External node " + i + " has no taxon");
            }

            final String trait = getTraitString(tree, node, taxon, traitName);
            writer.println("    <taxon id=\"" + escapeXml(taxon.getId()) + "\">");
            writer.println("        <attr name=\"" + escapeXml(traitName) + "\">" +
                    escapeXml(trait) + "</attr>");
            writer.println("    </taxon>");
        }
        writer.println("</taxa>");
        writer.flush();
    }

    private static String getTraitString(final Tree tree,
                                         final NodeRef node,
                                         final Taxon taxon,
                                         final String traitName) throws XMLParseException {
        final Object taxonTrait = taxon.getAttribute(traitName);
        if (taxonTrait != null) {
            return taxonTrait.toString();
        }

        final Object nodeTrait = tree.getNodeAttribute(node, traitName);
        if (nodeTrait instanceof double[]) {
            return formatTrait((double[]) nodeTrait);
        }

        throw new XMLParseException("No simulated trait named '" + traitName +
                "' was found for taxon " + taxon.getId());
    }

    private static String formatTrait(final double[] values) {
        final StringBuilder builder = new StringBuilder();
        for (int i = 0; i < values.length; ++i) {
            if (i > 0) {
                builder.append(' ');
            }
            builder.append(values[i]);
        }
        return builder.toString();
    }

    private static String escapeXml(final String value) {
        final StringBuilder builder = new StringBuilder(value.length());
        for (int i = 0; i < value.length(); ++i) {
            final char c = value.charAt(i);
            switch (c) {
                case '&':
                    builder.append("&amp;");
                    break;
                case '<':
                    builder.append("&lt;");
                    break;
                case '>':
                    builder.append("&gt;");
                    break;
                case '"':
                    builder.append("&quot;");
                    break;
                case '\'':
                    builder.append("&apos;");
                    break;
                default:
                    builder.append(c);
                    break;
            }
        }
        return builder.toString();
    }

    public static final XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        @Override
        public String getParserName() {
            return GENERAL_OU_TRAIT_SIMULATOR;
        }

        @Override
        public Object parseXMLObject(final XMLObject xo) throws XMLParseException {
            final Tree simulatedTree = (Tree) ((AbstractXMLObjectParser) GeneralOUTreeSimulator.PARSER).parseXMLObject(xo);
            final String traitName = xo.getStringAttribute(GeneralOUTreeSimulator.TRAIT_NAME);
            final String taxaId = xo.getAttribute(TAXA_ID, "simulatedTaxa");

            final PrintWriter writer = XMLParser.getFilePrintWriter(xo, getParserName(), FILE_NAME);
            writeTaxaBlock(simulatedTree, traitName, taxaId, writer);
            writer.close();

            return simulatedTree;
        }

        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            return RULES;
        }

        @Override
        public String getParserDescription() {
            return "Simulates exact multivariate OU tip traits and writes them as a BEAST taxa block.";
        }

        @Override
        public Class getReturnType() {
            return Tree.class;
        }
    };

    private static final XMLSyntaxRule[] RULES = makeRules();

    private static XMLSyntaxRule[] makeRules() {
        final XMLSyntaxRule[] baseRules = GeneralOUTreeSimulator.PARSER.getSyntaxRules();
        final XMLSyntaxRule[] rules = new XMLSyntaxRule[baseRules.length + 3];
        rules[0] = new StringAttributeRule(FILE_NAME,
                "The file to write the simulated BEAST taxa block to.");
        rules[1] = AttributeRule.newBooleanRule(LoggerParser.ALLOW_OVERWRITE_LOG, true);
        rules[2] = new StringAttributeRule(TAXA_ID,
                "The id to use for the written taxa block.", true);
        System.arraycopy(baseRules, 0, rules, 3, baseRules.length);
        return rules;
    }
}
