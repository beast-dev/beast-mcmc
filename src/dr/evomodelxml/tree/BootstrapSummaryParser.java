/*
 * BootstrapSummaryParser.java
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

package dr.evomodelxml.tree;

import dr.evolution.io.Importer;
import dr.evomodel.tree.BootstrapSummary;
import dr.util.FileHelpers;
import dr.xml.*;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;

/**
 * Summarizes a bootstrap once its replicates have run, when parsed (BootstrapSummary):
 * <pre>
 * &lt;bootstrapSummary level="0.95" parameters="run.params.tsv" bootstrapParameters="run.bootstrap.params.tsv"
 *                   summary="run.summary.txt"&gt;
 *     &lt;treeFiles timeTree="run.timetree.nex" bootstrapTrees="run.bootstrap.trees" dates="run.dates.tsv"
 *                rates="run.rates.tsv"/&gt;
 * &lt;/bootstrapSummary&gt;
 * </pre>
 *
 * @author Xiang Ji
 * @author Marc A. Suchard
 */
public class BootstrapSummaryParser extends AbstractXMLObjectParser {

    public static final String BOOTSTRAP_SUMMARY = "bootstrapSummary";
    public static final String LEVEL = "level";
    public static final String PARAMETERS = "parameters";
    public static final String BOOTSTRAP_PARAMETERS = "bootstrapParameters";
    public static final String SUMMARY = "summary";
    public static final String TREE_FILES = "treeFiles";
    public static final String TIME_TREE = "timeTree";
    public static final String BOOTSTRAP_TREES = "bootstrapTrees";
    public static final String DATES = "dates";
    public static final String RATES = "rates";

    public String getParserName() {
        return BOOTSTRAP_SUMMARY;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        double level = xo.getAttribute(LEVEL, 0.95);
        if (!(level > 0.0 && level < 1.0)) {
            throw new XMLParseException("The level of " + BOOTSTRAP_SUMMARY + " must lie between 0 and 1");
        }
        List<BootstrapSummary.TreeFiles> trees = new ArrayList<BootstrapSummary.TreeFiles>();
        for (int i = 0; i < xo.getChildCount(); ++i) {
            if (xo.getChild(i) instanceof XMLObject && ((XMLObject) xo.getChild(i)).getName().equals(TREE_FILES)) {
                XMLObject child = (XMLObject) xo.getChild(i);
                trees.add(new BootstrapSummary.TreeFiles(file(child, TIME_TREE), file(child, BOOTSTRAP_TREES),
                        file(child, DATES), child.hasAttribute(RATES) ? file(child, RATES) : null));
            }
        }
        BootstrapSummary summary = new BootstrapSummary(file(xo, PARAMETERS), file(xo, BOOTSTRAP_PARAMETERS),
                xo.hasAttribute(SUMMARY) ? file(xo, SUMMARY) : null, trees, level);
        try {
            int replicates = summary.summarize();
            Logger.getLogger("dr.evomodel").info("bootstrapSummary: " + replicates + " replicates summarized");
        } catch (IOException e) {
            throw new XMLParseException("bootstrapSummary: " + e.getMessage());
        } catch (Importer.ImportException e) {
            throw new XMLParseException("bootstrapSummary cannot read a tree: " + e.getMessage());
        }
        return summary;
    }

    private static File file(XMLObject xo, String attribute) throws XMLParseException {
        return FileHelpers.getFile(xo.getStringAttribute(attribute));
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "Writes percentile intervals from bootstrap replicates into the original analysis's tables and time tree.";
    }

    public Class getReturnType() {
        return BootstrapSummary.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newDoubleRule(LEVEL, true),
            AttributeRule.newStringRule(PARAMETERS),
            AttributeRule.newStringRule(BOOTSTRAP_PARAMETERS),
            AttributeRule.newStringRule(SUMMARY, true),
            new ElementRule(TREE_FILES, new XMLSyntaxRule[]{
                    AttributeRule.newStringRule(TIME_TREE),
                    AttributeRule.newStringRule(BOOTSTRAP_TREES),
                    AttributeRule.newStringRule(DATES),
                    AttributeRule.newStringRule(RATES, true),
            }, 0, Integer.MAX_VALUE),
    };
}
