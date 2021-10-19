package dr.evomodel.treedatalikelihood.discrete;

import dr.evolution.alignment.PatternList;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.TreeTrait;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.inference.model.Statistic;
import dr.inference.operators.hmc.MassPreconditioner;
import dr.xml.Reportable;

/**
 * A statistic that computes the maximum likelihood estimates between sequences based on a SubstitutionModel CTMC
 *
 * @author Andy Magee
 * @author Marc A. Suchard
 */
public class SequenceDistanceStatistic extends Statistic.Abstract implements Reportable {

    public enum DistanceType {
        MAXIMIZED_DISTANCE("distance", "distanceTo"),
        LOG_LIKELIHOOD("likelihood", "lnL");

        DistanceType(String name, String label) {
            this.name = name;
            this.label = label;
        }

        public String getName() {
            return name;
        }

        public String getLabel() { return label; }

        private String name;
        private String label;
    }

    public SequenceDistanceStatistic(TreeDataLikelihood treeLike, SubstitutionModel subsModel, PatternList patterns,
                                     boolean treeSeqAncestral, DistanceType type) {
        this.treeDataLikelihood = treeLike;
        this.substitutionModel = subsModel;
        this.patternList = patterns;
        this.treeSequenceIsAncestral = treeSeqAncestral;
        this.type = type;
    }
//    public void setTree(Tree tree) {
//        this.tree = tree;
//    }
//
//    public Tree getTree() {
//        return tree;
//    }

    public int getDimension() {
        return patternList.getTaxonCount();
    }

    public String getDimensionName(int i) {
        StringBuffer sb = new StringBuffer();
        sb.append(type.getLabel());
        sb.append("(");
        sb.append(patternList.getTaxonId(i));
        sb.append(")");
        return sb.toString();
    }

    public String getStatisticName() {
        return NAME;
    }

    /**
     * @return the statistic
     */
    public double getStatisticValue(int dim) {
        // Eventually we may want to enable this for other nodes
        NodeRef node = treeDataLikelihood.getTree().getRoot();

        return 42;
    }

    @Override
    public String getReport() {

        NodeRef node = treeDataLikelihood.getTree().getRoot();

        StringBuilder sb = new StringBuilder("sequenceDistanceStatistic Report\n\n");

        TreeTrait[] traits = treeDataLikelihood.getTreeTraits();
        for (TreeTrait trait : traits) {
            System.err.println(trait.getTraitName().toString());
        }

        System.err.println("treeDataLikelihood.getTreeTraits().length = " + treeDataLikelihood.getTreeTraits().length);

        sb.append("42");

        sb.append("\n\n");

        sb.append(patternList.getPatternWeights().length);

        sb.append("\n\n");

        sb.append(treeDataLikelihood.getTreeTraits().length);

        sb.append("\n\n");

        sb.append(treeDataLikelihood.getDataLikelihoodDelegate().getModelCount());

        sb.append("\n\n");

        double[] mat = new double[substitutionModel.getFrequencyModel().getFrequencyCount()*substitutionModel.getFrequencyModel().getFrequencyCount()];
        substitutionModel.getTransitionProbabilities(1, mat);
        sb.append(mat);

        treeDataLikelihood.getTreeTraits();
        sb.append("\n\n");

        return sb.toString();
    }

    private TreeDataLikelihood treeDataLikelihood = null;
    private PatternList patternList = null;
    private SubstitutionModel substitutionModel = null;
    boolean treeSequenceIsAncestral;
    private DistanceType type;

}
