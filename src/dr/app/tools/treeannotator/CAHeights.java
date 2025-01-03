/*
 * CAHeights.java
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

package dr.app.tools.treeannotator;

import dr.evolution.io.Importer;
import dr.evolution.tree.MutableTree;

import java.io.IOException;
import java.util.BitSet;

public class CAHeights {
    private final TreeAnnotator treeAnnotator;

    public CAHeights(TreeAnnotator treeAnnotator) {
        this.treeAnnotator = treeAnnotator;
    }

    // very inefficient, but Java wonderful bitset has no subset op
    // perhaps using bit iterator would be faster, I can't br bothered.
    static boolean isSubSet(BitSet x, BitSet y) {
        y = (BitSet) y.clone();
        y.and(x);
        return y.equals(x);
    }

    boolean setTreeHeightsByCA(MutableTree targetTree, final String inputFileName, final int burnin)
            throws IOException, Importer.ImportException {

//        long startTime = System.currentTimeMillis();
//
//        TreeAnnotator.progressStream.println("Setting node heights...");
//        TreeAnnotator.progressStream.println("0              25             50             75            100");
//        TreeAnnotator.progressStream.println("|--------------|--------------|--------------|--------------|");
//
//        final FileReader fileReader = new FileReader(inputFileName);
//        final NexusImporter importer = new NexusImporter(fileReader, true);
//
//        // this call increments the clade counts and it shouldn't
//        // this is remedied with removeClades call after while loop below
//        CladeSystem cladeSystem = new CladeSystem(treeAnnotator, targetTree);
//        final int nClades = cladeSystem.getCladeMap().size();
//
//        // allocate posterior tree nodes order once
//        int[] postOrderList = new int[nClades];
//        BitSet[] ctarget = new BitSet[nClades];
//        BitSet[] ctree = new BitSet[nClades];
//
//        for (int k = 0; k < nClades; ++k) {
//            ctarget[k] = new BitSet();
//            ctree[k] = new BitSet();
//        }
//
//        cladeSystem.getTreeCladeCodes(targetTree, ctarget);
//
//        // temp collecting heights inside loop allocated once
//        double[] hs = new double[nClades];
//
//        // heights total sum from posterior trees
//        double[] ths = new double[nClades];
//
//        treeAnnotator.setTotalTreesUsed(0);
//
//        int counter = 0;
//        while (importer.hasTree()) {
//            final Tree tree = importer.importNextTree();
//
//            if (counter >= burnin) {
//                TreeUtils.preOrderTraversalList(tree, postOrderList);
//                cladeSystem.getTreeCladeCodes(tree, ctree);
//                for (int k = 0; k < nClades; ++k) {
//                    int j = postOrderList[k];
//                    for (int i = 0; i < nClades; ++i) {
//                        if (isSubSet(ctarget[i], ctree[j])) {
//                            hs[i] = tree.getNodeHeight(tree.getNode(j));
//                        }
//                    }
//                }
//                for (int k = 0; k < nClades; ++k) {
//                    ths[k] += hs[k];
//                }
//                treeAnnotator.setTotalTreesUsed(treeAnnotator.getTotalTreesUsed() + 1);
//            }
//            if (counter > 0 && counter % reportStepSize == 0) {
//                TreeAnnotator.progressStream.print("*");
//                TreeAnnotator.progressStream.flush();
//            }
//            counter++;
//
//        }
//        cladeSystem.removeClades(targetTree, targetTree.getRoot(), true);
//        for (int k = 0; k < nClades; ++k) {
//            ths[k] /= treeAnnotator.getTotalTreesUsed();
//            final NodeRef node = targetTree.getNode(k);
//            targetTree.setNodeHeight(node, ths[k]);
//        }
//        fileReader.close();
//
//        long timeElapsed = (System.currentTimeMillis() - startTime) / 1000;
//        TreeAnnotator.progressStream.println("* [" + timeElapsed + " secs]");
//        TreeAnnotator.progressStream.println();

        return true;
    }
}