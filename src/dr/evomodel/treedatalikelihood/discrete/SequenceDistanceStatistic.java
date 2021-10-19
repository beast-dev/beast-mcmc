package dr.evomodel.treedatalikelihood.discrete;

import dr.evolution.alignment.PatternList;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evolution.tree.TreeUtils;
import dr.evolution.util.TaxonList;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treelikelihood.AncestralStateBeagleTreeLikelihood;
import dr.inference.model.Statistic;
import dr.inference.operators.hmc.MassPreconditioner;
import dr.math.UnivariateFunction;
import dr.math.UnivariateMinimum;
import dr.xml.Reportable;

import java.util.Set;

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
                                     boolean treeSeqAncestral,
                                     TaxonList mrcaTaxa, DistanceType type) throws TreeUtils.MissingTaxonException {
        this.asrLikelihood = asrLike;
        this.substitutionModel = subsModel;
        this.patternList = patterns;
        this.treeSequenceIsAncestral = treeSeqAncestral;
        this.type = type;
        this.tree = asrLikelihood.getTreeModel();
        this.leafSet = (mrcaTaxa != null) ? TreeUtils.getLeavesForTaxa(tree, mrcaTaxa) : null;
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
        NodeRef node = (leafSet != null) ?  TreeUtils.getCommonAncestorNode(tree, leafSet) : tree.getRoot();
//
//        int[] rootState = asrLikelihood.getStatesForNode(tree, node);
//
//        for (int s : rootState) {
//            System.err.println(s);
//        }
//
//        //asrLikelihood.getPatternsList().getTaxonIndex(taxa[dim]) should work when/if we have taxon list in our input

        return optimizeBranchLength(dim);
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

    private double optimizeBranchLength(int taxonIndex) {

        // Eventually we may want to enable this for other node
        int[] rootState = asrLikelihood.getStatesForNode(asrLikelihood.getTreeModel(),asrLikelihood.getTreeModel().getRoot());

        //asrLikelihood.getPatternsList().getTaxonIndex(taxa[dim]) should work when/if we have taxon list in our input

        int nStates = substitutionModel.getFrequencyModel().getFrequencyCount();

        UnivariateFunction f = new UnivariateFunction() {
            @Override
            public double evaluate(double argument) {

                double[] tpmFlat = new double[nStates*nStates];
                substitutionModel.getTransitionProbabilities(argument, tpmFlat);

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
                return -lnL;
            }

            @Override
            public double getLowerBound() {
                return 0;
            }

            @Override
            public double getUpperBound() {
                // TODO: should use some constant times the tree length in substitutions
                return 10.0;
            }
        };

        UnivariateMinimum minimum = new UnivariateMinimum();

        double x = minimum.findMinimum(f);

        double val;

        if ( type.getName() == "likelihood" ) {
            val = -minimum.fminx;
        } else {
            val = minimum.minx;
        }
        System.err.println("Used " + minimum.numFun + " evaluations to find minimum at " + minimum.minx + " with function value " + minimum.fminx + " and curvature " + minimum.f2minx);
        return val;
    }

    private AncestralStateBeagleTreeLikelihood asrLikelihood = null;
    private PatternList patternList = null;
    private SubstitutionModel substitutionModel = null;
    boolean treeSequenceIsAncestral;
    private final DistanceType type;
    private final Set<String> leafSet;
    private final Tree tree;

}
