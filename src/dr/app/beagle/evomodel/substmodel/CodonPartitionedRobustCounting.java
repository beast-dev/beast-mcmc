package dr.app.beagle.evomodel.substmodel;

import dr.app.beagle.evomodel.treelikelihood.AncestralStateBeagleTreeLikelihood;
import dr.evolution.tree.Tree;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.BranchAttributeProvider;
//import dr.evolution.datatype.DataType;
//import dr.evolution.datatype.Nucleotides;
//import dr.app.beagle.evomodel.substmodel.CodonLabeling;
import dr.inference.loggers.Loggable;
import dr.inference.loggers.LogColumn;
import dr.inference.loggers.NumberColumn;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Variable;
import dr.inference.model.Parameter;

import java.util.List;
import java.util.ArrayList;

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

public class CodonPartitionedRobustCounting extends AbstractModel implements BranchAttributeProvider, Loggable {

    public CodonPartitionedRobustCounting(String name, Tree tree,
                                          AncestralStateBeagleTreeLikelihood partition0,
                                          AncestralStateBeagleTreeLikelihood partition1,
                                          AncestralStateBeagleTreeLikelihood partition2,
                                          CodonLabeling codonLabeling) {

        super(name);
        this.tree = tree;

        this.partition0 = partition0;
        this.partition1 = partition1;
        this.partition2 = partition2;
        this.codonLabeling = codonLabeling;

        substModel0 = partition0.getSubstitutionModel();
        substModel1 = partition1.getSubstitutionModel();
        substModel2 = partition2.getSubstitutionModel();

        List<SubstitutionModel> substModels = new ArrayList<SubstitutionModel>(3);
//        List<DataType> dataTypes = new ArrayList<DataType>(3);

        substModels.add(substModel0);
//        dataTypes.add(Nucleotides.INSTANCE);
        substModels.add(substModel1);
//        dataTypes.add(Nucleotides.INSTANCE);
        substModels.add(substModel2);
//        dataTypes.add(Nucleotides.INSTANCE);

        productChain = new ProductChainSubstitutionModel("codonLabeling", //dataTypes,
                substModels);

        addModel(substModel0);
        addModel(substModel1);
        addModel(substModel2);
    }

    public int getCountsForBranch(NodeRef child) {

        // Get child node reconstructed sequence
        final int[] childSeq0 = partition0.getStatesForNode(tree, child);
        final int[] childSeq1 = partition1.getStatesForNode(tree, child);
        final int[] childSeq2 = partition2.getStatesForNode(tree, child);

        // Get parent node reconstructed sequence
        final NodeRef parent = tree.getParent(child);
        final int[] parentSeq0 = partition0.getStatesForNode(tree, parent);
        final int[] parentSeq1 = partition1.getStatesForNode(tree, parent);
        final int[] parentSeq2 = partition2.getStatesForNode(tree, parent);

        final int numCodons = childSeq0.length;
        final int[] childCodon  = new int[3];
        final int[] parentCodon = new int[3];

        int count = 0;

        for (int i = 0; i < numCodons; i++) {

            // Construct this child and parent codon
            childCodon[0] = childSeq0[i];
            childCodon[1] = childSeq1[i];
            childCodon[2] = childSeq2[i];

            parentCodon[0] = parentSeq0[i];
            parentCodon[1] = parentSeq1[i];
            parentCodon[2] = parentSeq2[i];

            if (DEBUG) {

                System.err.println("Computing robust counts for " +
                    parentCodon[0] + parentCodon[1] + parentCodon[2] + " -> " +
                    childCodon[0] + childCodon[1] + childCodon[2]);
            }

            // TODO Compute robust counts
        }

        return count;
    }

    public String getBranchAttributeLabel() {
        return codonLabeling.getText();
    }

    public String getAttributeForBranch(Tree tree, NodeRef node) {
        return String.valueOf(getCountsForBranch(node));
    }

    public int getRobustCount() {
        return 0;
    }


    public LogColumn[] getColumns() {
        return new LogColumn[]{
                new RobustCountColumn(codonLabeling.getText())
        };
    }

    protected void handleModelChangedEvent(Model model, Object object, int index) {
        // Do nothing
    }

    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        // Do nothing
    }

    protected void storeState() {
        // Do nothing
    }

    protected void restoreState() {
        // Do nothing
    }

    protected void acceptState() {
        // Do nothing
    }

    private class RobustCountColumn extends NumberColumn {

        public RobustCountColumn(String label) {
            super(label);
        }

        public double getDoubleValue() {
            return getRobustCount();
        }
    }

    private final AncestralStateBeagleTreeLikelihood partition0;
    private final AncestralStateBeagleTreeLikelihood partition1;
    private final AncestralStateBeagleTreeLikelihood partition2;

    private final SubstitutionModel substModel0;
    private final SubstitutionModel substModel1;
    private final SubstitutionModel substModel2;

    private final ProductChainSubstitutionModel productChain;

    private final CodonLabeling codonLabeling;
    private final Tree tree;

    private static final boolean DEBUG = true;
}
