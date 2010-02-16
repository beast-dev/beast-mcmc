package dr.app.beagle.evomodel.substmodel;

import dr.app.beagle.evomodel.treelikelihood.AncestralStateBeagleTreeLikelihood;
import dr.evolution.tree.Tree;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.BranchAttributeProvider;
import dr.inference.markovjumps.CodonLabeling;

/**
 * @author Marc A. Suchard
 * @author Vladimir Minin
 *         <p/>
 *         A class for implementing robust counting for synonymous and nonsynonymous changes in BEAST using BEAGLE
 *         This work is supported by NSF grant 0856099
 *         <p/>
 *         O'Brien JD, Minin VN and Suchard MA (2009) Learning to count: robust estimates for labeled distances between
 *         molecular sequences. Molecular Biology and Evolution, 26, 801-814 
 */

public class CodonPartitionedRobustCounting implements BranchAttributeProvider {

    public static final String SYN_COUNTS = "countS";
    public static final String NONSYN_COUNTS = "countN";

    public CodonPartitionedRobustCounting(Tree tree,
                                          AncestralStateBeagleTreeLikelihood partition1,
                                          AncestralStateBeagleTreeLikelihood partition2,
                                          AncestralStateBeagleTreeLikelihood partition3,
                                          CodonLabeling codonLabeling) {
        this.tree = tree;
        this.partition1 = partition1;
        this.partition2 = partition2;
        this.partition3 = partition3;
        this.codonLabeling = codonLabeling;
    }

    public int getCountsForBranch(NodeRef child) {

        final int[] childSeq0 = partition1.getStatesForNode(tree, child);
        final int[] childSeq1 = partition2.getStatesForNode(tree, child);
        final int[] childSeq2 = partition3.getStatesForNode(tree, child);

        final NodeRef parent = tree.getParent(child);
        final int[] parentSeq0 = partition1.getStatesForNode(tree, parent);
        final int[] parentSeq1 = partition2.getStatesForNode(tree, parent);
        final int[] parentSeq2 = partition3.getStatesForNode(tree, parent);

        final int numCodons = childSeq0.length;
        final int[] childCodon  = new int[3];
        final int[] parentCodon = new int[3];

        int count = 0;

        for (int i = 0; i < numCodons; i++) {
            childCodon[0] = childSeq0[i];
            childCodon[1] = childSeq1[i];
            childCodon[2] = childSeq2[i];

            parentCodon[0] = parentSeq0[i];
            parentCodon[1] = parentSeq1[i];
            parentCodon[2] = parentSeq2[i];

        }

        return count;
    }

    public String getBranchAttributeLabel() {
        return SYN_COUNTS;
    }

    public String getAttributeForBranch(Tree tree, NodeRef node) {
        return String.valueOf(getCountsForBranch(node));
    }

//    public String[] getNodeAttributeLabel() {
//        return new String[] { SYN_COUNTS, NONSYN_COUNTS };
//    }
//
//    public String[] getAttributeForNode(Tree tree, NodeRef node) {
//
//        if (this.tree != tree) {
//            throw new RuntimeException("Can only compute codon partition counts for a codon partitioned tree");
//        }
//
//        return new String[] {
//                String.valueOf(getCountsForBranch(node)),
//                String.valueOf(getCountsForBranch(node))
//        };
//    }

    private final AncestralStateBeagleTreeLikelihood partition1;
    private final AncestralStateBeagleTreeLikelihood partition2;
    private final AncestralStateBeagleTreeLikelihood partition3;
    private final CodonLabeling codonLabeling;
    private final Tree tree;

}
