/*
 * TestParsimony.java
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

package dr.evolution.parsimony;

import dr.evolution.alignment.*;
import dr.evolution.io.NexusImporter;
import dr.evolution.tree.*;
import dr.evolution.datatype.PairedDataType;

import java.io.FileReader;

/**
 * @author rambaut
 *         Date: Jun 21, 2005
 *         Time: 11:41:32 PM
 */
public class TestParsimony {

    public final static String alignmentFileName = "ABLV.G.nuc";
    public final static String treeFileName = "ABLV.G.tree";

    public static void main(String[] args) throws Exception {

        // read alignments
        NexusImporter importer = new NexusImporter(new FileReader(alignmentFileName));
        Alignment alignment = importer.importAlignment();
       // alignment = new GapStrippedAlignment(alignment);

        // read trees
        Tree tree = null;
        importer = new NexusImporter(new FileReader(treeFileName));
        try {
            tree = importer.importTree(alignment);
            //((FlexibleTree)tree).resolveTree();

        } catch (Exception e) {
            System.err.println("Exception attempting to read tree: " + e.getMessage());
            System.exit(0);
        }

        //testParsimony(new FitchParsimony(alignment, false), alignment, tree);
        System.out.println("single:");
        testParsimony(new SankoffParsimony(alignment), alignment, tree);
        System.out.println("paired:");
        testPairedSites(alignment, tree);
    }

    static void testPairedSites(SiteList siteList, Tree tree) {

        ParsimonyCriterion singleParsimony = new SankoffParsimony(siteList);
        double[] singleScores = singleParsimony.getSiteScores(tree);

        PairedDataType pairedDataType = new PairedDataType(siteList.getDataType());
        for (int i = 0; i < siteList.getSiteCount(); i++) {
            if (singleScores[i] > 0) {
                SimpleSiteList pairedSites = new SimpleSiteList(pairedDataType, siteList);

                int[] siteIndices = new int[siteList.getSiteCount()];
                int u = 0;

                for (int j = i + 1; j < siteList.getSiteCount(); j++) {
                    if (singleScores[j] > 0) {

                        int[] pairedPattern = new int[siteList.getTaxonCount()];

                        int[] pattern1 = siteList.getSitePattern(i);
                        int[] pattern2 = siteList.getSitePattern(j);

                        for (int k = 0; k < pairedPattern.length; k ++) {
                            pairedPattern[k] = pairedDataType.getState(pattern1[k], pattern2[k]);
                        }
                        pairedSites.addPattern(pairedPattern);

                        siteIndices[u] = j;
                        u++;
                    }
                }

                //testParsimony(new FitchParsimony(pairedSites, false), pairedSites, tree);
                ParsimonyCriterion parsimony = new SankoffParsimony(pairedSites);
                double[] scores = parsimony.getSiteScores(tree);

                for (int k = 0; k < pairedSites.getSiteCount(); k++) {
                    if (scores[k] > 1 /*&& isEpistatic(parsimony, pairedDataType, k, tree, tree.getRoot())*/) {
                        System.out.println("site {" + Integer.toString(i + 1) + ", " + Integer.toString(siteIndices[k] + 1)  + "}: " +  scores[k] + " steps");
                        showChanges(parsimony, k, tree, tree.getRoot());
                    }
                }
            }
        }

    }

    static void testParsimony(ParsimonyCriterion parsimony, SiteList siteList, Tree tree) {
        double[] scores = parsimony.getSiteScores(tree);

        for (int i = 0; i < siteList.getSiteCount(); i++) {
            if (scores[i] > 1) {
                System.out.println("site" + i  + " (" +  scores[i] + " steps):");
                showChanges(parsimony, i, tree, tree.getRoot());
            }
        }

    }

    static void showChanges(ParsimonyCriterion parsimony, int site, Tree tree, NodeRef node) {

        if (node != tree.getRoot()) {
            int stateB = parsimony.getStates(tree, node)[site];
            int stateA = parsimony.getStates(tree, tree.getParent(node))[site];
            if (stateA != stateB) {
                if (tree.isExternal(node)) {
                    System.out.println("\tnode(" + tree.getParent(node).getNumber() + " -> " + tree.getNodeTaxon(node).getId() + "):\tstate" + "(" +  stateA + " -> " + stateB + ")");
                } else {
                    System.out.println("\tnode(" + tree.getParent(node).getNumber() + " -> " + node.getNumber() + "):\tstate" + "(" +  stateA + " -> " + stateB + ")");
                }
            }
        }
        for (int j = 0; j < tree.getChildCount(node); j++) {
            showChanges(parsimony, site, tree, tree.getChild(node, j));
        }

    }

}
