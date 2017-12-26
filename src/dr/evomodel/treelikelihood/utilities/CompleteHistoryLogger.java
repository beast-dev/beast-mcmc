/*
 * CompleteHistoryLogger.java
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

package dr.evomodel.treelikelihood.utilities;

import dr.evomodel.treelikelihood.MarkovJumpsBeagleTreeLikelihood;
import dr.evomodel.treelikelihood.MarkovJumpsTraitProvider;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.inference.loggers.LogColumn;
import dr.inference.loggers.Loggable;
import dr.inference.markovjumps.StateHistory;
import dr.util.Citable;
import dr.util.Citation;
import dr.util.CommonCitations;

import java.awt.*;
import java.io.Serializable;
import java.util.*;
import java.util.List;
import java.util.logging.Logger;

/**
 * A class to conveniently log a complete state history of a continuous-time Markov chain along a tree
 * simulated using the Uniformization Method
 * <p/>
 * This work is supported by NSF grant 0856099
 * <p/>
 * Minin VN and Suchard MA (2008) Counting labeled transitions in continous-time Markov models of evolution.
 * Journal of Mathematical Biology, 56, 391-412.
 * <p/>
 * Rodrigue N, Philippe H and Lartillot N (2006) Uniformization for sampling realizations of Markov processes:
 * applications to Bayesian implementations of codon substitution models. Bioinformatics, 24, 56-62.
 * <p/>
 * Hobolth A and Stone E (2009) Simulation from endpoint-conditioned, continuous-time Markov chains on a finite
 * state space, with applications to molecular evolution. Annals of Applied Statistics, 3, 1204-1231.
 *
 * @author Marc A. Suchard
 * @author Philippe Lemey
 * @author Andrew Rambaut
 */
public class CompleteHistoryLogger implements Loggable, Citable {

    public static final String TOTAL_COUNT_NAME = "totalChangeCount";
    public static final String COMPLETE_HISTORY_NAME = "completeHistory";

    public CompleteHistoryLogger(MarkovJumpsTraitProvider treeLikelihood, HistoryFilter filter, boolean internal, boolean external) {
        this.tree = treeLikelihood.getTreeModel();
        this.patternCount = treeLikelihood.getPatternCount();
        this.internal = internal;
        this.external = external;

        treeTraitHistory = new TreeTrait[patternCount];
        for (int site = 0; site < patternCount; ++site) {
            String traitName = (patternCount == 1) ? MarkovJumpsBeagleTreeLikelihood.HISTORY : MarkovJumpsBeagleTreeLikelihood.HISTORY + "_" + (site + 1);
            treeTraitHistory[site] = treeLikelihood.getTreeTrait(traitName);
            if (treeTraitHistory[site] == null) {
                throw new RuntimeException("Tree '" + treeLikelihood.getId() + "' does not have a complete history trait at site " + (site + 1));
            }
        }

        treeTraitCount = treeLikelihood.getTreeTrait(MarkovJumpsBeagleTreeLikelihood.TOTAL_COUNTS);

        if (treeTraitCount == null) {
            throw new RuntimeException("No sum");
        }

        if (filter == null) {
            this.filter = new HistoryFilter.Default();
        } else {
            this.filter = filter;
            Logger.getLogger("dr.app.beagle").info("\tWith filter: " + filter.getDescription() + "\n");
        }


    }

    private static int parseListString(String listString, int currentOffset, List list) {

        while (currentOffset < listString.length()) {
            // Skip leading separators
            if (listString.startsWith(",", currentOffset)) {
                currentOffset++;
            }

            if (listString.startsWith("{", currentOffset)) {
                // Need to make a new list
                List newList = new ArrayList();
                currentOffset = parseListString(listString, currentOffset + 1, newList);
                list.add(newList);
            } else if (listString.startsWith("}", currentOffset)) {
                // the list is ended
                return currentOffset + 1;
            } else {
                // Parse until next ',' or '}'
                int nextComma = listString.indexOf(",", currentOffset);
                int nextClose = listString.indexOf("}", currentOffset);
                if (nextComma < 0) nextComma = listString.length() - 1;
                if (nextClose < 0) nextClose = listString.length() - 1;
                int nextOffset = Math.min(nextComma, nextClose);

                list.add(listString.substring(currentOffset, nextOffset).trim());
                currentOffset = nextOffset;
            }
        }
        return currentOffset;
    }

    public static Serializable parseValue(String value) {
        List nestedParse = new ArrayList<Object>();
        parseListString(value, 0, nestedParse);
        return parseValueObject(nestedParse.get(0));
    }

    private static Serializable parseValueObject(Object objValue) {
        if (objValue instanceof List) {
            List newObjList = (List) objValue;
            Object[] newObjects = new Object[newObjList.size()];
            for (int i = 0; i < newObjList.size(); ++i) {
                newObjects[i] = parseValueObject(newObjList.get(i));
            }
            return newObjects;
        }

        String value = (String) objValue;
        if (value.startsWith("#")) {
            // I am not sure whether this is a good idea but
            // I am going to assume that a # denotes an RGB colour
            try {
                return Color.decode(value.substring(1));
            } catch (NumberFormatException nfe1) {
                // not a colour
            }
        }

        if (value.equalsIgnoreCase("TRUE") || value.equalsIgnoreCase("FALSE")) {
            return Boolean.valueOf(value);
        }

        // Attempt to format the value as an integer
        try {
            return new Integer(value);
        } catch (NumberFormatException nfe1) {
            // not an integer
        }

        // Attempt to format the value as a double
        try {
            return new Double(value);
        } catch (NumberFormatException nfe2) {
            // not a double
        }

        // return the trimmed string
        return value;

    }


    public void setFilter(HistoryFilter filter) {
        this.filter = filter;
    }

    public LogColumn[] getColumns() {

        LogColumn[] columns = new LogColumn[1 + patternCount];
        columns[0] = new LogColumn.Abstract(TOTAL_COUNT_NAME) {

            @Override
            protected String getFormattedValue() {
                return treeTraitCount.getTraitString(tree, null);
            }
        };

        for (int site = 0; site < patternCount; ++site) {

            String name = (patternCount == 0) ? COMPLETE_HISTORY_NAME : COMPLETE_HISTORY_NAME + "_" + (site + 1);
            final int anonSite = site;
            columns[1 + site] = new LogColumn.Abstract(name) {

                @Override
                protected String getFormattedValue() {
                    boolean empty = true;
                    StringBuilder bf = new StringBuilder("{");
                    int count = 0;
                    for (int i = 0; i < tree.getNodeCount(); ++i) {
                        NodeRef node = tree.getNode(i);
                        if (!tree.isRoot(node)) {
                            if((tree.isExternal(node) && external) || (!tree.isExternal(node) && internal)){
                                NodeRef parent = tree.getParent(node);
                                double parentTime = tree.getNodeHeight(parent);
                                double childTime = tree.getNodeHeight(node);
                                double minTime = Math.min(parentTime, childTime);
                                double maxTime = Math.max(parentTime, childTime);
                                String trait = treeTraitHistory[anonSite].getTraitString(tree, node);
                                if (trait != null && trait.compareTo("{}") != 0) {
                                    Object[] changes = (Object[]) parseValue(trait);
                                    for (int j = 0; j < changes.length; ++j) {

                                        Object[] change = (Object[]) changes[j];
                                        int offset = (change.length == 4) ? 1 : 0;
                                        String source = (String) change[1 + offset];
                                        String dest = (String) change[2 + offset];
                                        double thisTime = (Double) change[0 + offset];
                                        if (thisTime < 0.0) {
                                            throw new RuntimeException("negative time");
                                        }
                                        if (thisTime > maxTime || thisTime < minTime) {
                                            throw new RuntimeException("Invalid simulation time");
                                        }

                                        // TODO Delegate to Filter(source, dest, thisTime).  If filtered then
                                        boolean filtered = filter.filter(source, dest, thisTime);
                                        if (filtered) {
                                            if (!empty) {
                                                bf.append(",");
                                            }
                                            StateHistory.addEventToStringBuilder(bf, source, dest,
                                                    thisTime, anonSite + 1);
                                            count++;
                                            empty = false;
                                        } else {
                                            // Do nothing
                                        }
                                    }
                                }
                            }
                        }
                    }
                    bf.append("}").append(" ").append(count);
                    return bf.toString();
                }
            };
        }
        return columns;
    }

    @Override
    public Citation.Category getCategory() {
        return Citation.Category.COUNTING_PROCESSES;
    }

    @Override
    public String getDescription() {
        return "Complete history logger";
    }

    public List<Citation> getCitations() {
        return Arrays.asList(//CommonCitations.LEMEY_2012, // TODO Find published Lemey paper
                CommonCitations.MININ_2008_FAST, CommonCitations.BLOOM_2013_STABILITY);
    }

    final private Tree tree;
    final private TreeTrait[] treeTraitHistory;
    final private TreeTrait treeTraitCount;
    final private int patternCount;
    final private boolean internal;
    final private boolean external;

    private HistoryFilter filter;
}
