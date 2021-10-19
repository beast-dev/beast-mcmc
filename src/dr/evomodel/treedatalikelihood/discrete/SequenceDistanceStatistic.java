package dr.evomodel.treedatalikelihood.discrete;

import dr.evolution.alignment.PatternList;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.TreeTrait;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treelikelihood.AncestralStateBeagleTreeLikelihood;
import dr.inference.model.Statistic;
import dr.xml.Reportable;

/**
 * A statistic that tracks the time of MRCA of a set of taxa
 *
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @version $Id: TMRCAStatistic.java,v 1.21 2005/07/11 14:06:25 rambaut Exp $
 */
public class SequenceDistanceStatistic extends Statistic.Abstract implements Reportable {

    public SequenceDistanceStatistic(AncestralStateBeagleTreeLikelihood asrLike, SubstitutionModel subsModel, PatternList patterns, boolean treeSeqAncestral, boolean reportDists) {
        this.asrLikelihood = asrLike;
        this.substitutionModel = subsModel;
        this.patternList = patterns;
        treeSequenceIsAncestral = treeSeqAncestral;
        reportDistances = reportDists;
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
        String dimname = "";
        if ( reportDistances ) {
            dimname += "distanceTo(";
        } else {
            dimname += "lnL(";
        }
        dimname += patternList.getTaxon(i);
        dimname += ")";
        return dimname;
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

        StringBuilder sb = new StringBuilder("sequenceDistanceStatistic Report\n\n");

//        TreeTrait[] traits = asrLikelihood.getTreeTraits();
//        for (TreeTrait trait : traits) {
//            System.err.println(trait.getTraitName().toString());
//        }
//
//        System.err.println("treeDataLikelihood.getTreeTraits().length = " + asrLikelihood.getTreeTraits().length);

        sb.append(getStatisticValue(1));

        sb.append("\n\n");

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
    boolean reportDistances;

}
