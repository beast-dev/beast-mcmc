/*
 * TreeLoggerParser.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
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
import dr.inference.loggers.Columns;
import dr.inference.loggers.LogFormatter;
import dr.inference.loggers.Loggable;
import dr.inference.loggers.TabDelimitedFormatter;
import dr.inference.model.Likelihood;
import dr.inference.model.Model;
import dr.inferencexml.loggers.LoggerParser;
import dr.xml.*;

import java.io.PrintWriter;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;

/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @author Karthik Gangavarapu
 */
public class TreeReportParser extends LoggerParser {

    public static final String REPORT_TREE = "reportTree";
    public static final String BRANCH_LENGTHS = "branchLengths";
    public static final String SUBSTITUTIONS = "substitutions";
    public static final String SORT_TRANSLATION_TABLE = "sortTranslationTable";
    public static final String DECIMAL_PLACES = "dp";
    public static final String ALLOW_OVERWRITE_LOG = "overwrite";

    public static final String FILTER_TRAITS = "traitFilter";
    public static final String TREE_TRAIT = "trait";
    public static final String NAME = "name";
    public static final String TAG = "tag";

    public String getParserName() {
        return REPORT_TREE;
    }

    protected String parseXMLParameters(XMLObject xo) throws XMLParseException
    {
        // reset this every time...
        branchRates = null;

        tree = (Tree) xo.getChild(Tree.class);

        sortTranslationTable = xo.getAttribute(SORT_TRANSLATION_TABLE, true);

        boolean substitutions = xo.getAttribute(BRANCH_LENGTHS, "").equals(SUBSTITUTIONS);

        List<TreeTraitProvider> ttps = new ArrayList<TreeTraitProvider>();

        // ttps2 are for TTPs that are not specified within a Trait element. These are only
        // included if not already added through a trait element to avoid duplication of
        // (in particular) the BranchRates which is required for substitution trees.
        List<TreeTraitProvider> ttps2 = new ArrayList<TreeTraitProvider>();

        for (int i = 0; i < xo.getChildCount(); i++) {
            Object cxo = xo.getChild(i);

            // This needs to be refactored into using a TreeTrait if Colouring is resurrected...
//            if (cxo instanceof TreeColouringProvider) {
//                final TreeColouringProvider colouringProvider = (TreeColouringProvider) cxo;
//                baps.add(new BranchAttributeProvider() {
//
//                    public String getBranchAttributeLabel() {
//                        return "deme";
//                    }
//
//                    public String getAttributeForBranch(Tree tree, NodeRef node) {
//                        TreeColouring colouring = colouringProvider.getTreeColouring(tree);
//                        BranchColouring bcol = colouring.getBranchColouring(node);
//                        StringBuilder buffer = new StringBuilder();
//                        if (bcol != null) {
//                            buffer.append("{");
//                            buffer.append(bcol.getChildColour());
//                            for (int i = 1; i <= bcol.getNumEvents(); i++) {
//                                buffer.append(",");
//                                buffer.append(bcol.getBackwardTime(i));
//                                buffer.append(",");
//                                buffer.append(bcol.getBackwardColourAbove(i));
//                            }
//                            buffer.append("}");
//                        }
//                        return buffer.toString();
//                    }
//                });
//
//            } else

            if (cxo instanceof TreeTraitProvider) {
                if (xo.hasAttribute(FILTER_TRAITS)) {
                    String[] matches = ((String) xo.getAttribute(FILTER_TRAITS)).split("[\\s,]+");
                    TreeTraitProvider ttp = (TreeTraitProvider) cxo;
                    TreeTrait[] traits = ttp.getTreeTraits();
                    List<TreeTrait> filteredTraits = new ArrayList<TreeTrait>();
                    for (String match : matches) {
                        for (TreeTrait trait : traits) {
                            if (trait.getTraitName().startsWith(match)) {
                                filteredTraits.add(trait);
                            }
                        }
                    }
                    if (filteredTraits.size() > 0) {
                        ttps2.add(new TreeTraitProvider.Helper(filteredTraits));
                    }

                } else {
                    // Add all of them
                    ttps2.add((TreeTraitProvider) cxo);
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

                        ttps.add(new TreeTraitProvider.Helper(tag, new TreeTrait() {

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
                        TreeTrait[] traits = ttp.getTreeTraits();
                        List<TreeTrait> filteredTraits = new ArrayList<TreeTrait>();
                        for (String match : matches) {
                            for (TreeTrait trait : traits) {
                                if (trait.getTraitName().startsWith(match)) {
                                    filteredTraits.add(trait);
                                }
                            }
                        }
                        if (filteredTraits.size() > 0) {
                            ttps.add(new TreeTraitProvider.Helper(filteredTraits));
                        }

                    } else {
                        // neither named or filtered traits so just add them all
                        ttps.add(ttp);
                    }
                }
            }
            // Without this next block, branch rates get ignored :-(
            // BranchRateModels are now TreeTraitProviders so this is not needed.
//            if (cxo instanceof TreeTrait) {
//                final TreeTrait trait = (TreeTrait)cxo;
//                TreeTraitProvider ttp = new TreeTraitProvider() {
//                    public TreeTrait[] getTreeTraits() {
//                        return new TreeTrait[]  { trait };
//                    }
//
//                    public TreeTrait getTreeTrait(String key) {
//                        if (key.equals(trait.getTraitName())) {
//                            return trait;
//                        }
//                        return null;
//                    }
//                };
//                ttps.add(ttp);
//            }
            //}

        }

        // if we don't have any of the newer trait elements but we do have some tree trait providers
        // included directly then assume the user wanted to log these as tree traits (it may be an older
        // form XML).
//        if (ttps.size() == 0 && ttps2.size() > 0) {
//            ttps.addAll(ttps2);
//        }

        // The above code destroyed the logging of complete histories - which need to be logged by direct
        // inclusion of the codon partitioned robust counting TTP...
        if (ttps2.size() > 0) {
            ttps.addAll(ttps2);
        }

        if (substitutions) {
            branchRates = (BranchRates) xo.getChild(BranchRates.class);
        }
        if (substitutions && branchRates == null) {
            throw new XMLParseException("To log trees in units of substitutions a BranchRateModel must be provided");
        }

        // decimal places
        final int dp = xo.getAttribute(DECIMAL_PLACES, -1);
        if (dp != -1) {
            format = NumberFormat.getNumberInstance(Locale.ENGLISH);
            format.setMaximumFractionDigits(dp);
        }

        final PrintWriter pw = getLogFile(xo, getParserName());

        formatter = new TabDelimitedFormatter(pw);

        treeTraitProviders = ttps.toArray(new TreeTraitProvider[ttps.size()]);

        return TreeUtils.newick(tree, treeTraitProviders);
    }

    /**
     * @return an object based on the XML element it was passed.
     */
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        String report = parseXMLParameters(xo);

        formatter.logLine(report);

        return null;
    }

    protected Tree tree;
    protected String title;
    protected boolean nexusFormat;
    protected boolean sortTranslationTable;
    protected BranchRates branchRates;
    protected NumberFormat format = null;
    protected TreeLogger.LogUpon condition;
    protected LogFormatter formatter;
    protected TreeTraitProvider[] treeTraitProviders;

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
