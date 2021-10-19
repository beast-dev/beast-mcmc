package dr.evomodel.treedatalikelihood.discrete;

import dr.evolution.alignment.PatternList;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evolution.tree.TreeUtils;
import dr.evolution.util.TaxonList;
import dr.evomodel.branchratemodel.AbstractBranchRateModel;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.branchratemodel.StrictClockBranchRates;
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

        public double extractResultForType( double[] results ) {
            return name == "distance" ? results[0] : results[1];
        }

        private String name;
        private String label;
    }

    public SequenceDistanceStatistic(AncestralStateBeagleTreeLikelihood asrLike,
                                     SubstitutionModel subsModel,
                                     BranchRateModel branchRates,
                                     PatternList patterns,
                                     boolean treeSeqAncestral,
                                     TaxonList distTaxa,
                                     TaxonList mrcaTaxa,
                                     DistanceType type) throws TreeUtils.MissingTaxonException {
        this.asrLikelihood = asrLike;
        this.substitutionModel = subsModel;
        this.branchRates = branchRates;
        this.patternList = patterns;
        this.treeSequenceIsAncestral = treeSeqAncestral;
        this.type = type;
        this.tree = asrLikelihood.getTreeModel();
        this.distanceTaxa = mrcaTaxa;
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
        return type.getLabel() +
                "(" +
                patternList.getTaxonId(i) +
                ")";
    }

    public String getStatisticName() {
        return NAME;
    }

    /**
     * @return the statistic
     */
    public double getStatisticValue(int dim) {
        double[] optimized = optimizeBranchLength(dim);
        return type.extractResultForType(optimized);
    }

    @Override
    public String getReport() {

        NodeRef node = (leafSet != null) ?  TreeUtils.getCommonAncestorNode(tree, leafSet) : tree.getRoot();

        StringBuilder sb = new StringBuilder("sequenceDistanceStatistic Report\n\n");

        for (int i=0; i < patternList.getTaxonCount(); i++) {
            String source = treeSequenceIsAncestral ? "node " + node.getNumber() : "taxon" + patternList.getTaxonId(i);
            String target = treeSequenceIsAncestral ? "taxon " + patternList.getTaxonId(i) : "node " + node.getNumber();
            sb.append("distance (in calendar time) from " + source + " to " + target + " is " + getStatisticValue(i) + "\n");
        }
        sb.append("\n\n");

        return sb.toString();
    }

    private double[] optimizeBranchLength(int taxonIndex) {
        NodeRef node = (leafSet != null) ?  TreeUtils.getCommonAncestorNode(tree, leafSet) : tree.getRoot();
        int[] nodeState = asrLikelihood.getStatesForNode(tree,node);

        //asrLikelihood.getPatternsList().getTaxonIndex(taxa[dim]) should work when/if we have taxon list in our input

        int nStates = substitutionModel.getFrequencyModel().getFrequencyCount();

        // could consider getting from asrLikelihood, probably, at the cost of an additional taxon list but removing need for patterns argument
        // TODO: scale by (assumed fixed) branch rate (second-order todo: what branch rate to use if not fixed? throw error in mean time)
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
                if ( treeSequenceIsAncestral ) {
                    for (int i=0; i<nodeState.length; i++) {
                        from = nodeState[i];
                        to = patternList.getPatternState(taxonIndex,i);
                        lnL += tpm[from][to];
                    }
                } else {
                    for (int i=0; i<nodeState.length; i++) {
                        to = nodeState[i];
                        from = patternList.getPatternState(taxonIndex,i);
                        lnL += tpm[from][to];
                    }
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

        // MAS: should delegate via something like: val = type.getReturnValue(minimum);
        double results[] = {minimum.minx/branchRates.getBranchRate(tree,node),-minimum.fminx};

//        System.err.println("Used " + minimum.numFun + " evaluations to find minimum at " + minimum.minx + " with function value " + minimum.fminx + " and curvature " + minimum.f2minx);

        return results;
    }

    private AncestralStateBeagleTreeLikelihood asrLikelihood = null;
    private BranchRateModel branchRates = null;
    private PatternList patternList = null;
    private SubstitutionModel substitutionModel = null;
    boolean treeSequenceIsAncestral;
    private final TaxonList distanceTaxa;
    private final Set<String> leafSet;
    private final Tree tree;
    private final DistanceType type;
}
