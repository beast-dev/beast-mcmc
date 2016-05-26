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

import dr.evolution.colouring.TreeColouringProvider;
import dr.evolution.tree.*;
import dr.evomodel.tree.TreeLogger;
import dr.inference.loggers.LogFormatter;
import dr.inference.loggers.Loggable;
import dr.inference.loggers.TabDelimitedFormatter;
import dr.inference.model.Likelihood;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inferencexml.loggers.LoggerParser;
import dr.util.Identifiable;
import dr.xml.*;

import java.io.PrintWriter;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;

/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 */
public class TreeLoggerParser extends LoggerParser {

    public static final String LOG_TREE = "logTree";
    public static final String NEXUS_FORMAT = "nexusFormat";
    //    public static final String USING_RATES = "usingRates";
    public static final String BRANCH_LENGTHS = "branchLengths";
    public static final String TIME = "time";
    public static final String SUBSTITUTIONS = "substitutions";
    public static final String SORT_TRANSLATION_TABLE = "sortTranslationTable";
    public static final String MAP_NAMES = "mapNamesToNumbers";
    public static final String DECIMAL_PLACES = "dp";
    //    public static final String NORMALISE_MEAN_RATE_TO = "normaliseMeanRateTo";

    public static final String FILTER_TRAITS = "traitFilter";
    public static final String TREE_TRAIT = "trait";
    public static final String NAME = "name";
    public static final String TAG = "tag";

    public String getParserName() {
        return LOG_TREE;
    }

    protected void parseXMLParameters(XMLObject xo) throws XMLParseException
    {
        // reset this every time...
        branchRates = null;

        tree = (Tree) xo.getChild(Tree.class);

        title = xo.getAttribute(TITLE, "");

        nexusFormat = xo.getAttribute(NEXUS_FORMAT, false);

        sortTranslationTable = xo.getAttribute(SORT_TRANSLATION_TABLE, true);

        boolean substitutions = xo.getAttribute(BRANCH_LENGTHS, "").equals(SUBSTITUTIONS);

        List<TreeAttributeProvider> taps = new ArrayList<TreeAttributeProvider>();
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

            if (cxo instanceof Likelihood) {
                final Likelihood likelihood = (Likelihood) cxo;
                taps.add(new TreeAttributeProvider() {

                    public String[] getTreeAttributeLabel() {
                        return new String[] {"lnP"};
                    }

                    public String[] getAttributeForTree(Tree tree) {
                        return new String[] {Double.toString(likelihood.getLogLikelihood())};
                    }
                });
            }

            if (cxo instanceof TreeAttributeProvider) {
                taps.add((TreeAttributeProvider) cxo);
            }
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

            // be able to put arbitrary statistics in as tree attributes
            if (cxo instanceof Loggable) {
                final Loggable loggable = (Loggable) cxo;
                taps.add(new TreeAttributeProvider() {

                    public String[] getTreeAttributeLabel() {
                        String[] labels = new String[loggable.getColumns().length];
                        for (int i = 0; i < loggable.getColumns().length; i++) {
                            labels[i] = loggable.getColumns()[i].getLabel();
                        }
                        return labels;
                    }

                    public String[] getAttributeForTree(Tree tree) {
                        String[] values = new String[loggable.getColumns().length];
                        for (int i = 0; i < loggable.getColumns().length; i++) {
                            values[i] = loggable.getColumns()[i].getFormatted();
                        }
                        return values;
                    }
                });
            }

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

        // logEvery of zero only displays at the end
        logEvery = xo.getAttribute(LOG_EVERY, 0);

//        double normaliseMeanRateTo = xo.getAttribute(NORMALISE_MEAN_RATE_TO, Double.NaN);

        // decimal places
        final int dp = xo.getAttribute(DECIMAL_PLACES, -1);
        if (dp != -1) {
            format = NumberFormat.getNumberInstance(Locale.ENGLISH);
            format.setMaximumFractionDigits(dp);
        }

        final PrintWriter pw = getLogFile(xo, getParserName());

        formatter = new TabDelimitedFormatter(pw);

        treeAttributeProviders = new TreeAttributeProvider[taps.size()];
        taps.toArray(treeAttributeProviders);
        treeTraitProviders = new TreeTraitProvider[ttps.size()];
        ttps.toArray(treeTraitProviders);

        // I think the default should be to have names rather than numbers, thus the false default - AJD
        // I think the default should be numbers - using names results in larger files and end user never
        // sees the numbers anyway as any software loading the nexus files does the translation - JH
        mapNames = xo.getAttribute(MAP_NAMES, true);

        condition = logEvery == 0 ? (TreeLogger.LogUpon) xo.getChild(TreeLogger.LogUpon.class) : null;
    }

    /**
     * @return an object based on the XML element it was passed.
     */
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        parseXMLParameters(xo);

        TreeLogger logger = new TreeLogger(tree, branchRates,
                treeAttributeProviders, treeTraitProviders,
                formatter, logEvery, nexusFormat, sortTranslationTable, mapNames, format, condition/*,
                normaliseMeanRateTo*/);

        if (title != null) {
            logger.setTitle(title);
        }

        return logger;
    }

    protected Tree tree;
    protected String title;
    protected boolean nexusFormat;
    protected boolean sortTranslationTable;
    protected BranchRates branchRates;
    protected NumberFormat format = null;
    protected TreeLogger.LogUpon condition;
    protected boolean mapNames;
    protected LogFormatter formatter;
    protected TreeAttributeProvider[] treeAttributeProviders;
    protected TreeTraitProvider[] treeTraitProviders;
    protected int logEvery;

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newIntegerRule(LOG_EVERY, true),
            AttributeRule.newBooleanRule(ALLOW_OVERWRITE_LOG, true),
            new StringAttributeRule(FILE_NAME,
                    "The name of the file to send log output to. " +
                            "If no file name is specified then log is sent to standard output", true),
            new StringAttributeRule(TITLE, "The title of the log", true),
            AttributeRule.newBooleanRule(NEXUS_FORMAT, true,
                    "Whether to use the NEXUS format for the tree log"),
            AttributeRule.newBooleanRule(SORT_TRANSLATION_TABLE, true,
                    "Whether the translation table is sorted."),
            /*AttributeRule.newDoubleRule(NORMALISE_MEAN_RATE_TO, true,
                    "Value to normalise the mean rate to."),*/
            new StringAttributeRule(BRANCH_LENGTHS, "What units should the branch lengths be in",
                    new String[]{TIME, SUBSTITUTIONS}, true),
            AttributeRule.newStringRule(FILTER_TRAITS, true),
            AttributeRule.newBooleanRule(MAP_NAMES, true),
            AttributeRule.newIntegerRule(DECIMAL_PLACES, true),

            new ElementRule(Tree.class, "The tree which is to be logged"),
//            new ElementRule(BranchRates.class, true),
//            new ElementRule(TreeColouringProvider.class, true),
            new ElementRule(TREE_TRAIT,
                    new XMLSyntaxRule[] {
                            AttributeRule.newStringRule(NAME, false, "The name of the trait"),
                            AttributeRule.newStringRule(TAG, true, "The label of the trait to be used in the tree"),
                            new ElementRule(TreeAttributeProvider.class, "The trait provider")
                    }, 0, Integer.MAX_VALUE),
            new ElementRule(Likelihood.class, true),
            new ElementRule(Loggable.class, 0, Integer.MAX_VALUE),
            new ElementRule(TreeAttributeProvider.class, 0, Integer.MAX_VALUE),
            new ElementRule(TreeTraitProvider.class, 0, Integer.MAX_VALUE),
            new ElementRule(TreeLogger.LogUpon.class, true)
    };

    public String getParserDescription() {
        return "Logs a tree to a file";
    }

    public String getExample() {
        final String name = getParserName();
        return
                "<!-- The " + name + " element takes a treeModel to be logged -->\n" +
                        "<" + name + " " + LOG_EVERY + "=\"100\" " + FILE_NAME + "=\"log.trees\" "
                        + NEXUS_FORMAT + "=\"true\">\n" +
                        "	<treeModel idref=\"treeModel1\"/>\n" +
                        "</" + name + ">\n";
    }

    public Class getReturnType() {
        return TreeLogger.class;
    }
}
