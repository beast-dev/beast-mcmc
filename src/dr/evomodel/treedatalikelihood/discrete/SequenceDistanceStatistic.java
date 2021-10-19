package dr.evomodel.treedatalikelihood.discrete;

import dr.evolution.alignment.PatternList;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.TreeTrait;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treelikelihood.AncestralStateBeagleTreeLikelihood;
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

    public SequenceDistanceStatistic(AncestralStateBeagleTreeLikelihood asrLike, SubstitutionModel subsModel, PatternList patterns,
                                     boolean treeSeqAncestral, DistanceType type) {
        this.asrLikelihood = asrLike;
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
        // Eventually we may want to enable this for other node
        int[] rootState = asrLikelihood.getStatesForNode(asrLikelihood.getTreeModel(),asrLikelihood.getTreeModel().getRoot());
        // Eventually we may want to enable this for other nodes
        for (int s : rootState) {
            System.err.println(s);
        }

        //asrLikelihood.getPatternsList().getTaxonIndex(taxa[dim]) should work when/if we have taxon list in our input

        int nStates = substitutionModel.getFrequencyModel().getFrequencyCount();

        double[] tpmFlat = new double[nStates*nStates];
        substitutionModel.getTransitionProbabilities(1, tpmFlat);

        // Make indexing easier in likelihood computation
        // This is really ln(P) and not P
        double[][] tpm = new double[nStates][nStates];
        for (int i=0; i < nStates; i++) {
            for (int j=0; j < nStates; j++) {
                tpm[i][j] = Math.log(tpmFlat[i*nStates+j]);
            }
        }


        int from,to = -1;
        double lnL = 0;
        for (int i=0; i<rootState.length; i++) {
            from = rootState[i];
            to = asrLikelihood.getPatternsList().getPatternState(0,i);
            lnL += tpm[from][to];
        }

        return lnL;
    }

    @Override
    public String getReport() {

//        NodeRef node = asrLikelihood.getTree().getRoot();
//
        StringBuilder sb = new StringBuilder("sequenceDistanceStatistic Report\n\n");
//
//        TreeTrait[] traits = asrLikelihood.getTreeTraits();
//        for (TreeTrait trait : traits) {
//            System.err.println(trait.getTraitName().toString());
//        }
//
//        System.err.println("treeDataLikelihood.getTreeTraits().length = " + asrLikelihood.getTreeTraits().length);

        sb.append(getStatisticValue(0));

        sb.append("\n\n");

//        sb.append(patternList.getPatternWeights().length);
//
//        sb.append("\n\n");
//
//        sb.append(asrLikelihood.getTreeTraits().length);
//
//        sb.append("\n\n");
//
//        sb.append(asrLikelihood.getDataLikelihoodDelegate().getModelCount());
//
//        sb.append("\n\n");
//
//        double[] mat = new double[substitutionModel.getFrequencyModel().getFrequencyCount()*substitutionModel.getFrequencyModel().getFrequencyCount()];
//        substitutionModel.getTransitionProbabilities(1, mat);
//        sb.append(mat);
//
//        asrLikelihood.getTreeTraits();
//        sb.append("\n\n");

        return sb.toString();
    }

    private AncestralStateBeagleTreeLikelihood asrLikelihood = null;
    private PatternList patternList = null;
    private SubstitutionModel substitutionModel = null;
    boolean treeSequenceIsAncestral;
    private DistanceType type;

}
