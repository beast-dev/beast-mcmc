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

import dr.evolution.colouring.TreeColouring;
import dr.evolution.colouring.TreeColouringProvider;
import dr.evolution.tree.Tree;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.transmission.TransmissionHistoryModel;
import dr.inference.loggers.LogFormatter;
import dr.inference.loggers.MCLogger;
import dr.inference.loggers.MLLogger;
import dr.inference.loggers.TabDelimitedFormatter;
import dr.inference.model.Likelihood;
import dr.xml.*;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintWriter;

/**
 * A logger that logs tree and clade frequencies.
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 *
 * @version $Id: TreeLogger.java,v 1.25 2006/09/05 13:29:34 rambaut Exp $
 */
public class TreeLogger extends MCLogger {

    public static final String LOG_TREE = "logTree";
    public static final String NEXUS_FORMAT = "nexusFormat";
    public static final String USING_RATES = "usingRates";
    public static final String BRANCH_LENGTHS = "branchLengths";
    public static final String TIME = "time";
    public static final String SUBSTITUTIONS = "substitutions";

    private Tree tree;
    private BranchRateModel branchRateModel = null;
    private String rateLabel;

    private TreeColouringProvider treeColouringProvider = null;
    private String colouringLabel;

    private TransmissionHistoryModel transmissionHistoryModel = null;
    private String hostLabel;

    private Likelihood likelihood = null;
    private String likelihoodLabel;

    private boolean nexusFormat = false;
    public boolean usingRates = false;
    public boolean substitutions = false;

    /**
     * Constructor
     */
    public TreeLogger(Tree tree, BranchRateModel branchRateModel, String rateLabel,
                      TreeColouringProvider treeColouringProvider, String colouringLabel,
                      Likelihood likelihood, String likelihoodLabel,
                      LogFormatter formatter, int logEvery, boolean nexusFormat, boolean substitutions) {

        super(formatter, logEvery);

        this.nexusFormat = nexusFormat;
        this.branchRateModel = branchRateModel;
        this.rateLabel = rateLabel;

        this.treeColouringProvider = treeColouringProvider;
        this.colouringLabel = colouringLabel;

        this.likelihood = likelihood;
        this.likelihoodLabel = likelihoodLabel;

        if (branchRateModel != null) {
            this.substitutions = substitutions;
        }
        this.tree = tree;
    }

    public void startLogging() {

        if (nexusFormat) {
            int taxonCount = tree.getTaxonCount();
            logLine("#NEXUS");
            logLine("");
            logLine("Begin taxa;");
            logLine("\tDimensions ntax=" + taxonCount + ";");
            logLine("\tTaxlabels");
            for (int i = 0; i < taxonCount; i++) {
                logLine("\t\t" + tree.getTaxon(i).getId());
            }
            logLine("\t\t;");
            logLine("End;");
            logLine("");
            logLine("Begin trees;");

            // This is needed if the trees use numerical taxon labels
            logLine("\tTranslate");
            for (int i = 0; i < taxonCount; i++) {
                int k = i + 1;
                if (k < taxonCount) {
                    logLine("\t\t" + k + " " + tree.getTaxonId(i) + ",");
                } else {
                    logLine("\t\t" + k + " " + tree.getTaxonId(i));
                }
            }
            logLine("\t\t;");
        }
    }

    public void log(int state) {

        if (logEvery <= 0 || ((state % logEvery) == 0)) {
            StringBuffer buffer = new StringBuffer("tree STATE_");
            buffer.append(state);
            if (likelihood != null) {
                buffer.append(" [&");
                buffer.append(likelihoodLabel);
                buffer.append("=");
                buffer.append(likelihood.getLogLikelihood());
                buffer.append("]");
            }

            buffer.append(" = [&R] ");

            if (substitutions) {
                Tree.Utils.newick(tree, tree.getRoot(), false, Tree.Utils.LENGTHS_AS_SUBSTITUTIONS,
                        branchRateModel, null, null, null, buffer);
            } else {
                if (treeColouringProvider != null) {
                    TreeColouring colouring = treeColouringProvider.getTreeColouring(tree);

                    Tree.Utils.newick(tree, tree.getRoot(), false, Tree.Utils.LENGTHS_AS_TIME,
                            branchRateModel, rateLabel, colouring, colouringLabel, buffer);
                } else {
                    Tree.Utils.newick(tree, tree.getRoot(), false, Tree.Utils.LENGTHS_AS_TIME,
                            branchRateModel, rateLabel, null, null, buffer);
                }
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

        public String getParserName() { return LOG_TREE; }

        /**
         * @return an object based on the XML element it was passed.
         */
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            Tree tree = (Tree)xo.getChild(Tree.class);

            String fileName = null;
            String title = null;
            boolean nexusFormat = false;

            String colouringLabel = "demes";
            String rateLabel = "rate";
            String likelihoodLabel = "lnP";

            if (xo.hasAttribute(TITLE)) {
                title = xo.getStringAttribute(TITLE);
            }

            if (xo.hasAttribute(FILE_NAME)) {
                fileName = xo.getStringAttribute(FILE_NAME);
            }

            if (xo.hasAttribute(NEXUS_FORMAT)) {
                nexusFormat = xo.getBooleanAttribute(NEXUS_FORMAT);
            }

            boolean substitutions = false;
            if (xo.hasAttribute(BRANCH_LENGTHS)) {
                substitutions = xo.getStringAttribute(BRANCH_LENGTHS).equals(SUBSTITUTIONS);
            }

            BranchRateModel branchRateModel = (BranchRateModel)xo.getChild(BranchRateModel.class);

            TreeColouringProvider treeColouringProvider = (TreeColouringProvider)xo.getChild(TreeColouringProvider.class);

            Likelihood likelihood = (Likelihood)xo.getChild(Likelihood.class);

            // logEvery of zero only displays at the end
            int logEvery = 1;

            if (xo.hasAttribute(LOG_EVERY)) {
                logEvery = xo.getIntegerAttribute(LOG_EVERY);
            }

            PrintWriter pw = null;

            if (fileName != null) {

                try {
                    File file = new File(fileName);
                    String name = file.getName();
                    String parent = file.getParent();

                    if (!file.isAbsolute()) {
                        parent = System.getProperty("user.dir");
                    }

//					System.out.println("Writing log file to "+parent+System.getProperty("path.separator")+name);
                    pw = new PrintWriter(new FileOutputStream(new File(parent, name)));
                } catch (FileNotFoundException fnfe) {
                    throw new XMLParseException("File '" + fileName + "' can not be opened for " + getParserName() + " element.");
                }
            } else {
                pw = new PrintWriter(System.out);
            }

            LogFormatter formatter = new TabDelimitedFormatter(pw);

            TreeLogger logger = new TreeLogger(tree,
                    branchRateModel, rateLabel,
                    treeColouringProvider, colouringLabel,
                    likelihood, likelihoodLabel,
                    formatter, logEvery, nexusFormat, substitutions);

            if (title != null) {
                logger.setTitle(title);
            }

            return logger;
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************
        public XMLSyntaxRule[] getSyntaxRules() { return rules; }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
                AttributeRule.newIntegerRule(LOG_EVERY),
                new StringAttributeRule(FILE_NAME,
                        "The name of the file to send log output to. " +
                                "If no file name is specified then log is sent to standard output", true),
                new StringAttributeRule(TITLE, "The title of the log", true),
                AttributeRule.newBooleanRule(NEXUS_FORMAT, true,
                        "Whether to use the NEXUS format for the tree log"),
                new StringAttributeRule(BRANCH_LENGTHS, "What units should the branch lengths be in", new String[] { TIME, SUBSTITUTIONS }, true),
                new ElementRule(Tree.class, "The tree which is to be logged"),
                new ElementRule(BranchRateModel.class, true),
                new ElementRule(TreeColouringProvider.class, true),
                new ElementRule(Likelihood.class, true)
        };

        public String getParserDescription() {
            return "Logs a tree to a file";
        }

        public String getExample() {
            return
                    "<!-- The " + getParserName() + " element takes a treeModel to be logged -->\n" +
                            "<" + getParserName() + " " + LOG_EVERY + "=\"100\" " + FILE_NAME + "=\"log.trees\" " + NEXUS_FORMAT + "=\"true\">\n" +
                            "	<treeModel idref=\"treeModel1\"/>\n" +
                            "</" + getParserName() + ">\n";
        }

        public Class getReturnType() { return MLLogger.class; }
    };

}