package dr.app.tools.newtreeannotator;

import dr.evolution.io.Importer;
import dr.evolution.io.NexusImporter;
import dr.evolution.tree.MutableTree;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;

import java.io.FileReader;
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