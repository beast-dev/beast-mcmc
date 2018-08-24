/*
 * TreeLogger.java
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

package dr.evomodel.tree;

import dr.app.tools.NexusExporter;
import dr.evolution.tree.*;
import dr.inference.loggers.LogFormatter;
import dr.inference.loggers.MCLogger;
import dr.util.Keywordable;

import java.text.NumberFormat;
import java.util.*;

/**
 * A logger that logs tree and clade frequencies.
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: TreeLogger.java,v 1.25 2006/09/05 13:29:34 rambaut Exp $
 */
public class TreeLogger extends MCLogger {

    private Tree tree;
	private BranchRates branchRates = null;

    private TreeAttributeProvider[] treeAttributeProviders;
    private TreeTraitProvider[] treeTraitProviders;

    private boolean nexusFormat = false;
    public boolean usingRates = false;
    public boolean substitutions = false;
    private final Map<String, Integer> idMap = new HashMap<String, Integer>();
    private final List<String> taxaIds = new ArrayList<String>();
    private boolean mapNames = true;

    /*private double normaliseMeanRateTo = Double.NaN;
    boolean normaliseMeanRate = false;*/

    private NumberFormat format;
    private LogUpon condition = null;

    /**
     * Interface to indicate when to log a tree
     */
    public interface LogUpon {
        /**
         *
         * @param state
         * @return  True if log tree of this state.
         */
       boolean logNow(long state);
    }

    public TreeLogger(Tree tree, LogFormatter formatter, int logEvery, boolean nexusFormat,
                      boolean sortTranslationTable, boolean mapNames) {

        this(tree, null, null, null, formatter, logEvery, nexusFormat, sortTranslationTable, mapNames, null, null/*, Double.NaN*/);
    }

    public TreeLogger(Tree tree, LogFormatter formatter, int logEvery, boolean nexusFormat,
                      boolean sortTranslationTable, boolean mapNames, NumberFormat format) {

        this(tree, null, null, null, formatter, logEvery, nexusFormat, sortTranslationTable, mapNames, format, null/*, Double.NaN*/);
    }

    public TreeLogger(Tree tree, BranchRates branchRates,
                      TreeAttributeProvider[] treeAttributeProviders,
                      TreeTraitProvider[] treeTraitProviders,
                      LogFormatter formatter, int logEvery, boolean nexusFormat,
                      boolean sortTranslationTable, boolean mapNames, NumberFormat format,
                      TreeLogger.LogUpon condition) {

        super(formatter, logEvery, false);

        this.condition = condition;

        /*this.normaliseMeanRateTo = normaliseMeanRateTo;
        if(!Double.isNaN(normaliseMeanRateTo)) {
            normaliseMeanRate = true;
        }*/

        this.nexusFormat = nexusFormat;
        // if not NEXUS, can't map names
        this.mapNames = mapNames && nexusFormat;

        this.branchRates = branchRates;

        this.treeAttributeProviders = treeAttributeProviders;
        this.treeTraitProviders = treeTraitProviders;

        if (this.branchRates != null) {
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

        this.format = format;
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
                logLine("\t\t" + cleanTaxonName(taxaId));
                }

            logLine("\t\t;");
            logLine("End;");
            logLine("");
            logLine("Begin trees;");

            if (tree instanceof Keywordable) {
                List<String> keywords = ((Keywordable)tree).getKeywords();
                if (keywords.size() > 0) {
                    StringBuffer sb = new StringBuffer("[keywords:");
                    for (String keyword: keywords) {
                        sb.append(" ");
                        sb.append(keyword);
                    }
                    sb.append("]");
                    logLine(sb.toString());
                }
            }

            if (mapNames) {
                // This is needed if the trees use numerical taxon labels
                logLine("\tTranslate");
                int k = 1;
                for (String taxaId : taxaIds) {
                    if (k < taxonCount) {
                        logLine("\t\t" + k + " " + cleanTaxonName(taxaId) + ",");
                    } else {
                        logLine("\t\t" + k + " " + cleanTaxonName(taxaId));
                    }
                    k += 1;
                }
                logLine("\t\t;");
            }
        }
    }

    private String cleanTaxonName(String taxaId) {
        if (taxaId.matches(NexusExporter.SPECIAL_CHARACTERS_REGEX)) {
            if (taxaId.contains("\'")) {
                if (taxaId.contains("\"")) {
                    throw new RuntimeException("Illegal taxon name - contains both single and double quotes");
                }

                return "\"" + taxaId + "\"";
            }

            return "\'" + taxaId + "\'";
        }
        return taxaId;
    }

    public void log(long state) {

        /*if(normaliseMeanRate) {
            NormaliseMeanTreeRate.analyze(tree, normaliseMeanRateTo);
        }*/

        final boolean doIt = condition != null ? condition.logNow(state) :
                    (logEvery < 0 || ((state % logEvery) == 0));

        if ( doIt ) {
            StringBuffer buffer = new StringBuffer("tree STATE_");
            buffer.append(state);
            if (treeAttributeProviders != null) {
                boolean hasAttribute = false;
                for (TreeAttributeProvider tap : treeAttributeProviders) {
                    String[] attributeLabel = tap.getTreeAttributeLabel();
                    String[] attributeValue = tap.getAttributeForTree(tree);
                    for (int i = 0; i < attributeLabel.length; i++) {
                        if (!hasAttribute) {
                            buffer.append(" [&");
                            hasAttribute = true;
                        } else {
                            buffer.append(",");
                        }
                        buffer.append(attributeLabel[i]);
                        buffer.append("=");
                        buffer.append(attributeValue[i]);
                    }
                }
                if (hasAttribute) {
                    buffer.append("]");
                }
            }

            buffer.append(" = [&R] ");

            if (substitutions) {
                TreeUtils.newick(tree, tree.getRoot(), false, TreeUtils.BranchLengthType.LENGTHS_AS_SUBSTITUTIONS,
                        format, branchRates, treeTraitProviders, idMap, buffer);
            } else {
                //System.out.println(treeTraitProviders.length);
                TreeUtils.newick(tree, tree.getRoot(), !mapNames, TreeUtils.BranchLengthType.LENGTHS_AS_TIME,
                        format, null, treeTraitProviders, idMap, buffer);
            }

            buffer.append(";");
            logLine(buffer.toString());
        }
    }

    public void stopLogging() {
        logLine("End;");
        super.stopLogging();
    }

    public Tree getTree() {
		return tree;
	}

	public void setTree(Tree tree) {
		this.tree = tree;
	}

	public TreeAttributeProvider[] getTreeAttributeProviders() {
        return this.treeAttributeProviders;
    }

	public TreeTraitProvider[] getTreeTraitProviders() {
        return this.treeTraitProviders;
    }

    public BranchRates getBranchRates() {
        return this.branchRates;
    }

}