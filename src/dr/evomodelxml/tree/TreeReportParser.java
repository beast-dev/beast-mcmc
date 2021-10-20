/*
 * TreeLoggerParser.java
 *
 * Copyright (c) 2002-2021 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodelxml.tree;

import dr.evolution.tree.*;
import dr.evomodel.tree.TreeLogger;
import dr.inference.loggers.LogFormatter;
import dr.inference.loggers.TabDelimitedFormatter;
import dr.inference.model.Likelihood;
import dr.inference.model.Model;
import dr.inferencexml.loggers.LoggerParser;
import dr.xml.*;

import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;

/**
 * @author Karthik Gangavarapu
 * @author Marc A. Suchard
 */
public class TreeReportParser extends LoggerParser {

    private static final String REPORT_TREE = "reportTree";
    private static final String BRANCH_LENGTHS = "branchLengths";
    private static final String SUBSTITUTIONS = "substitutions";
    private static final String ALLOW_OVERWRITE_LOG = "overwrite";
    private static final String FILTER_TRAITS = "traitFilter";
    private static final String TREE_TRAIT = "trait";
    private static final String NAME = "name";
    private static final String TAG = "tag";

    public String getParserName() {
        return REPORT_TREE;
    }

    private void AddFilterTraits(TreeTraitProvider ttp,
                                 List<TreeTraitProvider> ttpList,
                                 String[] matches) {

        TreeTrait[] traits = ttp.getTreeTraits();
        List<TreeTrait> filteredTraits = new ArrayList<>();
        for (String match : matches) {
            for (TreeTrait trait : traits) {
                if (trait.getTraitName().startsWith(match)) {
                    filteredTraits.add(trait);
                }
            }
        }
        if (filteredTraits.size() > 0) {
            ttpList.add(new TreeTraitProvider.Helper(filteredTraits));
        }
    }

    private String printTreeWithAttributes(XMLObject xo) throws XMLParseException {
        // reset this every time...

        Tree tree = (Tree) xo.getChild(Tree.class);

        boolean substitutions = xo.getAttribute(BRANCH_LENGTHS, "").equals(SUBSTITUTIONS);

        List<TreeTraitProvider> ttpList = new ArrayList<>();

        // ttps2 are for TTPs that are not specified within a Trait element. These are only
        // included if not already added through a trait element to avoid duplication of
        // (in particular) the BranchRates which is required for substitution trees.
        List<TreeTraitProvider> ttpList2 = new ArrayList<>();

        for (int i = 0; i < xo.getChildCount(); i++) {

            Object cxo = xo.getChild(i);

            if (cxo instanceof TreeTraitProvider) {
                if (xo.hasAttribute(FILTER_TRAITS)) {
                    String[] matches = ((String) xo.getAttribute(FILTER_TRAITS)).split("[\\s,]+");
                    TreeTraitProvider ttp = (TreeTraitProvider) cxo;
                    AddFilterTraits(ttp, ttpList2, matches);
                } else {
                    // Add all of them
                    ttpList2.add((TreeTraitProvider) cxo);
                }
            }
            if (cxo instanceof XMLObject) {
                XMLObject xco = (XMLObject)cxo;
                if (xco.getName().equals(TREE_TRAIT)) {

                    TreeTraitProvider ttp = (TreeTraitProvider)xco.getChild(TreeTraitProvider.class);

                    if (xco.hasAttribute(NAME)) {
                        // a specific named trait is required (optionally with a tag to name it in the tree file)

                        String name = xco.getStringAttribute(NAME);
                        final TreeTrait trait = ttp.getTreeTrait(name);

                        if (trait == null) {
                            String childName = "TreeTraitProvider";

                            if (ttp instanceof Likelihood) {
                                childName = ((Likelihood)ttp).prettyName();
                            } else  if (ttp instanceof Model) {
                                childName = ((Model)ttp).getModelName();
                            }

                            throw new XMLParseException("Trait named, " + name + ", not found for " + childName);
                        }

                        final String tag;
                        if (xco.hasAttribute(TAG)) {
                            tag = xco.getStringAttribute(TAG);
                        } else {
                            tag = name;
                        }

                        ttpList.add(new TreeTraitProvider.Helper(tag, new TreeTrait() {

                            public String getTraitName() {
                                return tag;
                            }

                            public Intent getIntent() {
                                return trait.getIntent();
                            }

                            public Class getTraitClass() {
                                return trait.getTraitClass();
                            }

                            public Object getTrait(Tree tree, NodeRef node) {
                                return trait.getTrait(tree, node);
                            }

                            public String getTraitString(Tree tree, NodeRef node) {
                                return trait.getTraitString(tree, node);
                            }

                            public boolean getLoggable() {
                                return trait.getLoggable();
                            }
                        }));
                    } else if (xo.hasAttribute(FILTER_TRAITS)) {
                        // else a filter attribute is given to ask for all traits that starts with a specific
                        // string

                        String[] matches = ((String) xo.getAttribute(FILTER_TRAITS)).split("[\\s,]+");
                        AddFilterTraits(ttp, ttpList2, matches);
                    } else {
                        // neither named or filtered traits so just add them all
                        ttpList.add(ttp);
                    }
                }
            }
        }

        if (ttpList2.size() > 0) {
            ttpList.addAll(ttpList2);
        }

        BranchRates branchRates = null;
        if (substitutions) {
            branchRates = (BranchRates) xo.getChild(BranchRates.class);
        }
        if (substitutions && branchRates == null) {
            throw new XMLParseException("To log trees in units of substitutions a BranchRateModel must be provided");
        }

        TreeTraitProvider[] treeTraitProviders = ttpList.toArray(new TreeTraitProvider[0]);

        return TreeUtils.newick(tree, treeTraitProviders);
    }

    /**
     * @return an object based on the XML element it was passed.
     */
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        String tree = printTreeWithAttributes(xo);
        final PrintWriter pw = getLogFile(xo, getParserName());

        LogFormatter formatter = new TabDelimitedFormatter(pw);
        formatter.logLine(tree);

        return null;
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "Logs a tree to a file";
    }

    public String getExample() {
        final String name = getParserName();
        return
                "<!-- The " + name + " element takes a treeModel to be reported -->\n" +
                        "<" + name + " " + FILE_NAME + "=\"log.trees\" " +
                        "	<treeModel idref=\"treeModel1\"/>\n" +
                        "</" + name + ">\n";
    }

    public Class getReturnType() {
        return TreeLogger.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newBooleanRule(ALLOW_OVERWRITE_LOG, true),
            new StringAttributeRule(FILE_NAME,
                    "The name of the file to send log output to. " +
                            "If no file name is specified then log is sent to standard output", true)
    };
}
