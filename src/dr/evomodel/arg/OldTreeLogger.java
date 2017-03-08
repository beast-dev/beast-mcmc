/*
 * OldTreeLogger.java
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

package dr.evomodel.arg;

import dr.evolution.colouring.TreeColouring;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeUtils;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.coalescent.structure.ColourSamplerModel;
import dr.inference.loggers.LogFormatter;
import dr.inference.loggers.MCLogger;
import dr.inference.model.Likelihood;

/**
 * A logger that logs tree and clade frequencies.
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: OldTreeLogger.java,v 1.25 2006/09/05 13:29:34 rambaut Exp $
 */
public class OldTreeLogger extends MCLogger {

    public static final String LOG_TREE = "logTree";
    public static final String NEXUS_FORMAT = "nexusFormat";
    public static final String USING_RATES = "usingRates";
    public static final String BRANCH_LENGTHS = "branchLengths";
    public static final String TIME = "time";
    public static final String SUBSTITUTIONS = "substitutions";

    // The following were in MCLogger; where did they go?

    public static final String LOG = "log";
    public static final String ECHO = "echo";
    public static final String ECHO_EVERY = "echoEvery";
    public static final String TITLE = "title";
    public static final String FILE_NAME = "fileName";
    public static final String FORMAT = "format";
    public static final String TAB = "tab";
    public static final String HTML = "html";
    public static final String PRETTY = "pretty";
    public static final String LOG_EVERY = "logEvery";

    public static final String COLUMNS = "columns";
    public static final String COLUMN = "column";
    public static final String LABEL = "label";
    public static final String SIGNIFICANT_FIGURES = "sf";
    public static final String DECIMAL_PLACES = "dp";
    public static final String WIDTH = "width";


    private Tree tree;
    private BranchRateModel branchRateModel = null;
    private String rateLabel;

    private ColourSamplerModel colourSamplerModel = null;
    private String colouringLabel;

    private Likelihood likelihood = null;
    private String likelihoodLabel;

    private boolean nexusFormat = false;
    public boolean usingRates = false;
    public boolean substitutions = false;

    /**
     * Constructor
     */
    public OldTreeLogger(Tree tree, BranchRateModel branchRateModel, String rateLabel,
                         ColourSamplerModel colourSamplerModel, String colouringLabel,
                         Likelihood likelihood, String likelihoodLabel,
                         LogFormatter formatter, int logEvery, boolean nexusFormat, boolean substitutions) {

        super(formatter, logEvery, false);

        this.nexusFormat = nexusFormat;
        this.branchRateModel = branchRateModel;
        this.rateLabel = rateLabel;

        this.colourSamplerModel = colourSamplerModel;
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

            if (!useTaxonLabels()) {
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
    }

    public void log(long state) {

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

            buffer.append(additionalInfo());

            buffer.append(" = [&R] ");

            TreeColouring colouring = null;
            if (colourSamplerModel != null) {
                colouring = colourSamplerModel.getTreeColouring();
            }

            Tree printTree = getPrintTree();

            if (substitutions) {
                TreeUtils.newick(printTree, printTree.getRoot(), useTaxonLabels(), TreeUtils.BranchLengthType.LENGTHS_AS_SUBSTITUTIONS,
                        null, branchRateModel, null, null, buffer);
            } else {
                TreeUtils.newick(printTree, printTree.getRoot(), useTaxonLabels(), TreeUtils.BranchLengthType.LENGTHS_AS_TIME,
                        null, branchRateModel, null, null, buffer);
            }

            buffer.append(";");
            logLine(buffer.toString());
        }
    }

    protected String additionalInfo() {
        return "";
    }

    protected Tree getPrintTree() {
        return tree;
    }

    protected Tree getTree() {
        return tree;
    }


    protected boolean useTaxonLabels() {
        return false;
    }

    public void stopLogging() {

        logLine("End;");
        super.stopLogging();
    }

}