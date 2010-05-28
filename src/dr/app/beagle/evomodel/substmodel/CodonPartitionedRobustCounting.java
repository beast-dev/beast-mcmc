package dr.app.beagle.evomodel.substmodel;

import dr.app.beagle.evomodel.treelikelihood.AncestralStateBeagleTreeLikelihood;
import dr.app.beagle.evomodel.sitemodel.SiteRateModel;
import dr.evolution.tree.Tree;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.BranchAttributeProvider;
import dr.evolution.datatype.Codons;
import dr.inference.loggers.Loggable;
import dr.inference.loggers.LogColumn;
import dr.inference.loggers.NumberColumn;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Variable;
import dr.inference.model.Parameter;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.tree.TreeModel;

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

    public CodonPartitionedRobustCounting(String name, TreeModel tree,
                                          AncestralStateBeagleTreeLikelihood[] partition,
                                          Codons codons,
                                          CodonLabeling codonLabeling) {
        this(name, tree, partition, codons, codonLabeling,
                StratifiedTraitOutputFormat.SUM_OVER_SITES, StratifiedTraitOutputFormat.SUM_OVER_SITES);

    }

    public CodonPartitionedRobustCounting(String name, TreeModel tree,
                                          AncestralStateBeagleTreeLikelihood[] partition,
                                          Codons codons,
                                          CodonLabeling codonLabeling,
                                          StratifiedTraitOutputFormat branchFormat,
                                          StratifiedTraitOutputFormat logFormat) {
        super(name);
        this.tree = tree;
        addModel(tree);

        if (partition.length != 3) {
            throw new RuntimeException("CodonPartition models require 3 partitions");
        }

        this.partition = partition;
        this.codonLabeling = codonLabeling;
        branchRateModel = partition[0].getBranchRateModel();
        addModel(branchRateModel);

        List<SubstitutionModel> substModelsList = new ArrayList<SubstitutionModel>(3);
        List<SiteRateModel> siteRateModelsList = new ArrayList<SiteRateModel>(3);

        numCodons = partition[0].getPatternWeights().length;

        for (int i = 0; i < 3; i++) {
            substModelsList.add(partition[i].getSubstitutionModel());
            siteRateModelsList.add(partition[i].getSiteRateModel());
            if (partition[i].getPatternWeights().length != numCodons) {
                throw new RuntimeException("All sequence lengths must be equal in CodonPartitionedRobustCounting");
            }
        }

        ProductChainSubstitutionModel productChainModel =
                new ProductChainSubstitutionModel("codonLabeling", substModelsList, siteRateModelsList);
        addModel(productChainModel);

        markovJumps = new MarkovJumpsSubstitutionModel(productChainModel);
        double[] synRegMatrix = CodonLabeling.getRegisterMatrix(codonLabeling, codons, true);
        markovJumps.setRegistration(synRegMatrix);

        condMeanMatrix = new double[64 * 64];

        this.branchFormat = branchFormat;
        this.logFormat = logFormat;

        computedCounts = new double[tree.getNodeCount()][]; // TODO Temporary until there exists a helper class
    }

    public double[] getExpectedCountsForBranch(NodeRef child) { // TODO This function will implement TraitProvider
        if (!countsKnown) {
            computeAllExpectedCounts();
        }
        return computedCounts[child.getNumber()];
    }

    private void computeAllExpectedCounts() {
        for (int i = 0; i < tree.getNodeCount(); i++) {
            NodeRef child = tree.getNode(i);
            if (!tree.isRoot(child)) {
                computedCounts[child.getNumber()] = computeExpectedCountsForBranch(child);
            }
        }
        countsKnown = true;
    }

    private double[] computeExpectedCountsForBranch(NodeRef child) {

        // Get child node reconstructed sequence
        final int[] childSeq0 = partition[0].getStatesForNode(tree, child);
        final int[] childSeq1 = partition[1].getStatesForNode(tree, child);
        final int[] childSeq2 = partition[2].getStatesForNode(tree, child);

        // Get parent node reconstructed sequence
        final NodeRef parent = tree.getParent(child);
        final int[] parentSeq0 = partition[0].getStatesForNode(tree, parent);
        final int[] parentSeq1 = partition[1].getStatesForNode(tree, parent);
        final int[] parentSeq2 = partition[2].getStatesForNode(tree, parent);

        double branchRateTime = branchRateModel.getBranchRate(tree, child) * tree.getBranchLength(child);

        markovJumps.computeCondStatMarkovJumps(branchRateTime, condMeanMatrix);

        double[] count = new double[numCodons];

        for (int i = 0; i < numCodons; i++) {

            // Construct this child and parent codon

            final int childState = getCanonicalState(childSeq0[i], childSeq1[i], childSeq2[i]);
            final int parentState = getCanonicalState(parentSeq0[i], parentSeq1[i], parentSeq2[i]);

//            final int vChildState = getVladimirState(childSeq0[i], childSeq1[i], childSeq2[i]);
//            final int vParentState = getVladimirState(parentSeq0[i], parentSeq1[i], parentSeq2[i]);

            final double codonCount = condMeanMatrix[parentState * 64 + childState];

            if (DEBUG) {

                System.err.println("Computing robust counts for " +
                        parentSeq0[i] + parentSeq1[i] + parentSeq2[i] + " -> " +
                        childSeq0[i] + childSeq1[i] + childSeq2[i] + " : " +
//                        vParentState + " -> " + vChildState + " = " + codonCount);
                        parentState + " -> " + childState + " = " + codonCount);
                System.exit(-1);
            }

            count[i] = codonCount;
        }

        return count;
    }

    public String getBranchAttributeLabel() {
        return codonLabeling.getText();
    }

    public String getAttributeForBranch(Tree tree, NodeRef node) {
        double[] counts = getExpectedCountsForBranch(node);
        final String returnString;

        if (branchFormat == StratifiedTraitOutputFormat.SUM_OVER_SITES) {
            double total = 0;
            for (double x : counts) {
                total += x;
            }
            returnString = String.valueOf(total);
        } else if (branchFormat == StratifiedTraitOutputFormat.PER_SITE) {
            if (counts.length == 1) {
                return String.valueOf(counts[0]);
            }
            StringBuffer sb = new StringBuffer("{");
            for (int i = 0; i < counts.length - 1; i++)
                sb.append(Double.toString(counts[i])).append(",");
            sb.append(Double.toString(counts[counts.length - 1])).append("}");
            returnString = sb.toString();
        } else {
            throw new RuntimeException("Not yet implemented.");
        }
        return returnString;
    }

    public double[] getRobustCountPerSite() { // TODO Move to StratifiedTraitLogger
        double[] count = new double[numCodons];

        for (int i = 0; i < tree.getNodeCount(); i++) {
            NodeRef node = tree.getNode(i);
            if (!tree.isRoot(node)) {
                addToMatrix(count, getExpectedCountsForBranch(node));
            }
        }
        return count;
    }

    public double getRobustCount() { // TODO Move to StratifiedTraitLogger

        double[] count = getRobustCountPerSite();
        double total = 0;
        for (double x : count) {
            total += x;
        }
        return total;
    }

    public double getRobustCount(int site) { // TODO Move to StratifiedTraitLogger

        if (site < 0 || site >= numCodons) {
            throw new RuntimeException("Invalid site #");
        }
        double total = 0;
        for (int i = 0; i < tree.getNodeCount(); i++) {
            NodeRef node = tree.getNode(i);
            if (!tree.isRoot(node)) {
                total += getExpectedCountsForBranch(node)[site];
            }
        }
        return total;
    }

    private void addToMatrix(double[] total, double[] summant) { // TODO Move to StratifiedTraitLogger
        final int length = summant.length;
        for (int i = 0; i < length; i++) {
            total[i] += summant[i];
        }
    }

    private int getCanonicalState(int i, int j, int k) {
        return i * 16 + j * 4 + k;
    }

//    private int getVladimirState(int i, int j, int k) {
//        if (i == 1) i = 2;
//        else if (i == 2) i = 1;
//
//        if (j == 1) j = 2;
//        else if (j == 2) j = 1;
//
//        if (k == 1) k = 2;
//        else if (k == 2) k = 1;
//
//        return i * 16 + j * 4 + k + 1;
//    }

    public LogColumn[] getColumns() {

        if (logFormat == StratifiedTraitOutputFormat.SUM_OVER_SITES) {
            return new LogColumn[]{
                    new RobustCountColumn(codonLabeling.getText())
            };
        } else {
            throw new RuntimeException("Not yet implemented");
        }
    }

    public int getDimension() {
        return numCodons;
    }

    public double getUnconditionedTraitValue() {
        return markovJumps.getMarginalRate() * getExpectedTreeLength();
    }

    private double getExpectedTreeLength() {
        double expectedTreeLength = 0;
        for (int i = 0; i < tree.getNodeCount(); i++) {
            NodeRef node = tree.getNode(i);
            if (!tree.isRoot(node)) {
                expectedTreeLength += branchRateModel.getBranchRate(tree, node)
                        * tree.getBranchLength(node);
            }
        }
        return expectedTreeLength;
    }

    protected void handleModelChangedEvent(Model model, Object object, int index) {
        countsKnown = false;
    }

    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        countsKnown = false;
    }

    protected void storeState() {
        // Do nothing
    }

    protected void restoreState() {
        countsKnown = false;
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

    private final AncestralStateBeagleTreeLikelihood[] partition;
    private final MarkovJumpsSubstitutionModel markovJumps;
    private final BranchRateModel branchRateModel;

    private final CodonLabeling codonLabeling;
    private final Tree tree;

    private final StratifiedTraitOutputFormat branchFormat;
    private final StratifiedTraitOutputFormat logFormat;

    private final double[] condMeanMatrix;

    private int numCodons;

    private static final boolean DEBUG = false;

    private boolean countsKnown;
    private double[][] computedCounts; // TODO Temporary storage until generic TreeTraitProvider/Helpers are finished
}
