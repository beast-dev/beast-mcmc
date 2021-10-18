package dr.evomodel.treedatalikelihood.discrete;

import dr.evolution.alignment.PatternList;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evolution.tree.TreeUtils;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;
import dr.evomodel.siteratemodel.SiteRateModel;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evomodel.tree.TreeStatistic;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.inference.model.Statistic;
import dr.math.matrixAlgebra.Vector;
import dr.oldevomodel.sitemodel.SiteModel;
import dr.xml.Reportable;

import java.util.Arrays;
import java.util.Set;

/**
 * A statistic that tracks the time of MRCA of a set of taxa
 *
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @version $Id: TMRCAStatistic.java,v 1.21 2005/07/11 14:06:25 rambaut Exp $
 */
public class SequenceDistanceStatistic extends Statistic.Abstract implements Reportable {

    public SequenceDistanceStatistic(TreeDataLikelihood treeLike, SubstitutionModel subsModel, PatternList patterns, boolean treeSeqAncestral, boolean reportDists) {
        this.treeDataLikelihood = treeLike;
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
    boolean reportDistances;

}
