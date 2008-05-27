/*
 * TreeLogger.java
 *
 * Copyright (C) 2002-2006 Alexei Drummond and Andrew Rambaut
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

package dr.evomodel.tree;

import dr.app.tools.NexusExporter;
import dr.evolution.colouring.BranchColouring;
import dr.evolution.colouring.TreeColouring;
import dr.evolution.colouring.TreeColouringProvider;
import dr.evolution.tree.*;
import dr.inference.loggers.LogFormatter;
import dr.inference.loggers.MCLogger;
import dr.inference.loggers.MLLogger;
import dr.inference.loggers.TabDelimitedFormatter;
import dr.inference.model.Likelihood;
import dr.xml.*;

import java.io.PrintWriter;
import java.util.*;

/**
 * A logger that logs tree and clade frequencies.
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: TreeLogger.java,v 1.25 2006/09/05 13:29:34 rambaut Exp $
 */
public class TreeLogger extends MCLogger {

    public static final String LOG_TREE = "logTree";
    public static final String NEXUS_FORMAT = "nexusFormat";
    public static final String USING_RATES = "usingRates";
    public static final String BRANCH_LENGTHS = "branchLengths";
    public static final String TIME = "time";
    public static final String SUBSTITUTIONS = "substitutions";
    public static final String SORT_TRANSLATION_TABLE = "sortTranslationTable";
    public static final String MAP_NAMES = "mapNamesToNumbers";

    private Tree tree;
    private BranchRateController branchRateProvider = null;

    private TreeAttributeProvider[] treeAttributeProviders;
    private NodeAttributeProvider[] nodeAttributeProviders;
    private BranchAttributeProvider[] branchAttributeProviders;

    private boolean nexusFormat = false;
    public boolean usingRates = false;
    public boolean substitutions = false;
    private Map<String, Integer> idMap = new HashMap<String, Integer>();
    private List<String> taxaIds = new ArrayList<String>();
    private boolean mapNames = true;



    public TreeLogger(Tree tree, BranchRateController branchRateProvider,
                      TreeAttributeProvider[] treeAttributeProviders,
                      NodeAttributeProvider[] nodeAttributeProviders,
                      BranchAttributeProvider[] branchAttributeProviders,
                      LogFormatter formatter, int logEvery, boolean nexusFormat,
                      boolean sortTranslationTable, boolean mapNames) {

        super(formatter, logEvery);

        this.nexusFormat = nexusFormat;
        this.mapNames = mapNames;

        this.branchRateProvider = branchRateProvider;

        this.treeAttributeProviders = treeAttributeProviders;
        this.nodeAttributeProviders = nodeAttributeProviders;

        this.branchAttributeProviders = branchAttributeProviders;

        if (this.branchRateProvider != null) {
            this.substitutions = true;
        }
        this.tree = tree;

        for (int i = 0; i < tree.getTaxonCount(); i++) {
            taxaIds.add(tree.getTaxon(i).getId());
        }
        if (sortTranslationTable) {
            Collections.sort(taxaIds);
        }

        int k = 1;
        for (String taxaId : taxaIds) {
            idMap.put(taxaId, k);
            k += 1;
        }
    }

    public void startLogging() {

        if (nexusFormat) {
            int taxonCount = tree.getTaxonCount();
            logLine("#NEXUS");
            logLine("");
            logLine("Begin taxa;");
            logLine("\tDimensions ntax=" + taxonCount + ";");
            logLine("\tTaxlabels");

            for (String taxaId : taxaIds) {
                if (taxaId.matches(NexusExporter.SPECIAL_CHARACTERS_REGEX)) {
                    taxaId = "'" + taxaId + "'";
                }
                logLine("\t\t" + taxaId);
            }

            logLine("\t\t;");
            logLine("End;");
            logLine("");
            logLine("Begin trees;");

            // This is needed if the trees use numerical taxon labels
            logLine("\tTranslate");
            int k = 1;
            for (String taxaId : taxaIds) {
                if (taxaId.matches(NexusExporter.SPECIAL_CHARACTERS_REGEX)) {
                    taxaId = "'" + taxaId + "'";
                }
                if (k < taxonCount) {
                    logLine("\t\t" + k + " " + taxaId + ",");
                } else {
                    logLine("\t\t" + k + " " + taxaId);
                }
                k += 1;
            }
            logLine("\t\t;");
        }
    }

    public void log(int state) {

        if (logEvery <= 0 || ((state % logEvery) == 0)) {
            StringBuffer buffer = new StringBuffer("tree STATE_");
            buffer.append(state);
            if (treeAttributeProviders != null) {
                boolean hasAttribute = false;
                for (TreeAttributeProvider tap : treeAttributeProviders) {
                    if (!hasAttribute) {
                        buffer.append(" [&");
                        hasAttribute = true;
                    } else {
                        buffer.append(",");
                    }
                    buffer.append(tap.getTreeAttributeLabel());
                    buffer.append("=");
                    buffer.append(tap.getAttributeForTree(tree));

                }
                if (hasAttribute) {
                    buffer.append("]");
                }
            }

            buffer.append(" = [&R] ");

            if (substitutions) {
                Tree.Utils.newick(tree, tree.getRoot(), false, Tree.Utils.LENGTHS_AS_SUBSTITUTIONS,
                        branchRateProvider, nodeAttributeProviders, branchAttributeProviders, idMap, buffer);
            } else {
                Tree.Utils.newick(tree, tree.getRoot(), !mapNames, Tree.Utils.LENGTHS_AS_TIME,
                        null, nodeAttributeProviders, branchAttributeProviders, idMap, buffer);
            }

            buffer.append(";");
            logLine(buffer.toString());
        }
    }

    public void stopLogging() {

        logLine("End;");
        super.stopLogging();
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return LOG_TREE;
        }

        /**
         * @return an object based on the XML element it was passed.
         */
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            Tree tree = (Tree) xo.getChild(Tree.class);


            String title = null;
            boolean nexusFormat = false;
            boolean sortTranslationTable = true;

            if (xo.hasAttribute(TITLE)) {
                title = xo.getStringAttribute(TITLE);
            }

            if (xo.hasAttribute(NEXUS_FORMAT)) {
                nexusFormat = xo.getBooleanAttribute(NEXUS_FORMAT);
            }

            if (xo.hasAttribute(SORT_TRANSLATION_TABLE)) {
                sortTranslationTable = xo.getBooleanAttribute(SORT_TRANSLATION_TABLE);
            }

            boolean substitutions = false;
            if (xo.hasAttribute(BRANCH_LENGTHS)) {
                substitutions = xo.getStringAttribute(BRANCH_LENGTHS).equals(SUBSTITUTIONS);
            }

            List<TreeAttributeProvider> taps = new ArrayList<TreeAttributeProvider>();
            List<NodeAttributeProvider> naps = new ArrayList<NodeAttributeProvider>();
            List<BranchAttributeProvider> baps = new ArrayList<BranchAttributeProvider>();

            for (int i = 0; i < xo.getChildCount(); i++) {
                Object cxo = xo.getChild(i);
                if (cxo instanceof TreeColouringProvider) {
                    final TreeColouringProvider colouringProvider = (TreeColouringProvider) cxo;
                    baps.add(new BranchAttributeProvider() {

                        public String getBranchAttributeLabel() {
                            return "deme";
                        }

                        public String getAttributeForBranch(Tree tree, NodeRef node) {
                            TreeColouring colouring = colouringProvider.getTreeColouring(tree);
                            BranchColouring bcol = colouring.getBranchColouring(node);
                            StringBuilder buffer = new StringBuilder();
                            if (bcol != null) {
                                buffer.append("{");
                                buffer.append(bcol.getChildColour());
                                for (int i = 1; i <= bcol.getNumEvents(); i++) {
                                    buffer.append(",");
                                    buffer.append(bcol.getBackwardTime(i));
                                    buffer.append(",");
                                    buffer.append(bcol.getBackwardColourAbove(i));
                                }
                                buffer.append("}");
                            }
                            return buffer.toString();
                        }
                    });

                } else if (cxo instanceof Likelihood) {
                    final Likelihood likelihood = (Likelihood) cxo;
                    taps.add(new TreeAttributeProvider() {

                        public String getTreeAttributeLabel() {
                            return "lnP";
                        }

                        public String getAttributeForTree(Tree tree) {
                            return Double.toString(likelihood.getLogLikelihood());
                        }
                    });

                } //else {
                if (cxo instanceof TreeAttributeProvider) {
                    taps.add((TreeAttributeProvider) cxo);
                }
                if (cxo instanceof NodeAttributeProvider) {
                    naps.add((NodeAttributeProvider) cxo);
                }
                if (cxo instanceof BranchAttributeProvider) {
                    baps.add((BranchAttributeProvider) cxo);
                }
                //}
            }
            BranchRateController branchRateProvider = null;
            if (substitutions) {
                branchRateProvider = (BranchRateController) xo.getChild(BranchRateController.class);
            }
            if (substitutions && branchRateProvider == null) {
                throw new XMLParseException("To log trees in units of substitutions a BranchRateModel must be provided");
            }

            // logEvery of zero only displays at the end
            int logEvery = 1;

            if (xo.hasAttribute(LOG_EVERY)) {
                logEvery = xo.getIntegerAttribute(LOG_EVERY);
            }
            PrintWriter pw = getLogFile(xo, getParserName());

            LogFormatter formatter = new TabDelimitedFormatter(pw);

            TreeAttributeProvider[] treeAttributeProviders = new TreeAttributeProvider[taps.size()];
            taps.toArray(treeAttributeProviders);
            NodeAttributeProvider[] nodeAttributeProviders = new NodeAttributeProvider[naps.size()];
            naps.toArray(nodeAttributeProviders);
            BranchAttributeProvider[] branchAttributeProviders = new BranchAttributeProvider[baps.size()];
            baps.toArray(branchAttributeProviders);

            boolean mapNames = true;
            if( xo.hasAttribute(MAP_NAMES) && !xo.getBooleanAttribute(MAP_NAMES) )
                mapNames = false;

            TreeLogger logger =
                    new TreeLogger(tree, branchRateProvider,
                            treeAttributeProviders, nodeAttributeProviders, branchAttributeProviders,
                            formatter, logEvery, nexusFormat, sortTranslationTable, mapNames);

            if (title != null) {
                logger.setTitle(title);
            }

            return logger;
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************
        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                AttributeRule.newIntegerRule(LOG_EVERY),
                new StringAttributeRule(FILE_NAME,
                        "The name of the file to send log output to. " +
                                "If no file name is specified then log is sent to standard output", true),
                new StringAttributeRule(TITLE, "The title of the log", true),
                AttributeRule.newBooleanRule(NEXUS_FORMAT, true,
                        "Whether to use the NEXUS format for the tree log"),
                AttributeRule.newBooleanRule(SORT_TRANSLATION_TABLE, true,
                        "Whether the translation table is sorted."),
                new StringAttributeRule(BRANCH_LENGTHS, "What units should the branch lengths be in",
                        new String[]{TIME, SUBSTITUTIONS}, true),
                new ElementRule(Tree.class, "The tree which is to be logged"),
                new ElementRule(BranchRateController.class, true),
                new ElementRule(TreeColouringProvider.class, true),
                new ElementRule(Likelihood.class, true),
                new ElementRule(TreeAttributeProvider.class, 0, Integer.MAX_VALUE),
                new ElementRule(NodeAttributeProvider.class, 0, Integer.MAX_VALUE),
                new ElementRule(BranchAttributeProvider.class, 0, Integer.MAX_VALUE),
                AttributeRule.newBooleanRule(MAP_NAMES,true)
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
            return MLLogger.class;
        }
    };

}