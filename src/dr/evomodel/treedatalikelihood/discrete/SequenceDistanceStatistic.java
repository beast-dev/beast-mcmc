package dr.evomodel.treedatalikelihood.discrete;

import dr.evolution.alignment.PatternList;
import dr.evolution.datatype.DataType;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeUtils;
import dr.evolution.util.TaxonList;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evomodel.treelikelihood.AncestralStateBeagleTreeLikelihood;
import dr.inference.model.Statistic;
import dr.math.UnivariateFunction;
import dr.math.UnivariateMinimum;
import dr.math.matrixAlgebra.Vector;
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
        MAXIMIZED_DISTANCE("distance", "distanceFrom") {
            public double extractResultForType(double[] results) {
                return results[0];
            }
        },
        LOG_LIKELIHOOD("likelihood", "lnL") {
            public double extractResultForType(double[] results) {
                return results[1];
            }
        };

        DistanceType(String name, String label) {
            this.name = name;
            this.label = label;
        }

        public String getName() {
            return name;
        }

        public String getLabel() {
            return label;
        }

        public abstract double extractResultForType(double[] results);

        private String name;
        private String label;
    }

    public SequenceDistanceStatistic(AncestralStateBeagleTreeLikelihood asrLike,
                                     SubstitutionModel subsModel,
                                     BranchRateModel branchRates,
                                     PatternList patterns,
                                     TaxonList mrcaTaxa,
                                     DistanceType type) throws TreeUtils.MissingTaxonException {
        this.asrLikelihood = asrLike;
        this.substitutionModel = subsModel;
        this.branchRates = branchRates;
        this.patternList = patterns;
        this.dataType = patternList.getDataType();
        this.type = type;
        this.tree = asrLikelihood.getTreeModel();
        this.leafSet = (mrcaTaxa != null) ? TreeUtils.getLeavesForTaxa(tree, mrcaTaxa) : null;

    }

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

        NodeRef node = (leafSet != null) ? TreeUtils.getCommonAncestorNode(tree, leafSet) : tree.getRoot();

        StringBuilder sb = new StringBuilder("sequenceDistanceStatistic Report\n\n");

        sb.append("dimension names: ");

        int n = getDimension();

        double[] values = new double[n];
        for (int i = 0; i < n; i++) {
            sb.append(getDimensionName(i));
            if (i != n - 1) {
                sb.append(" ");
            }
            values[i] = getStatisticValue(i);
        }
        sb.append("\n\n");
        sb.append("values: ");
        sb.append(new Vector(values));
        sb.append("\n\n");

        return sb.toString();
    }

    private double computeLogLikelihood(double distance, int[] taxonStates, int[] nodeStates) {
        // could consider getting from asrLikelihood, probably, at the cost of an additional taxon list but removing need for patterns argument
        int nStates = dataType.getStateCount();

        double[] tpm = new double[nStates * nStates];
        substitutionModel.getTransitionProbabilities(distance, tpm);
        double[] logTpm = new double[nStates * nStates];
        for (int i = 0; i < nStates * nStates; i++) {
            logTpm[i] = Math.log(tpm[i]);
        }

        double[] pi = substitutionModel.getFrequencyModel().getFrequencies();

        double lnL = 0.0;
        for (int s = 0; s < taxonStates.length; s++) {
            double sum = 0.0;
            int taxonState = taxonStates[s];
            int nodeState = nodeStates[s];
            if ( taxonState < nStates ) {
                lnL += logTpm[taxonState * nStates + nodeState];
            } else {
                for (int i = 0; i < nStates; i++) {
                    sum += tpm[i * nStates + nodeState] * pi[i];
                }
                lnL += Math.log(sum);
            }
        }
        return lnL;
    }

    private double[] optimizeBranchLength(int taxonIndex) {
        NodeRef node = (leafSet != null) ? TreeUtils.getCommonAncestorNode(tree, leafSet) : tree.getRoot();

        int[] nodeStates = asrLikelihood.getStatesForNode(tree,node);
        int[] taxonStates = new int[nodeStates.length];

        for (int i = 0; i < nodeStates.length; i++) {
            taxonStates[i] = patternList.getPatternState(taxonIndex,i);
        }

        UnivariateFunction f = new UnivariateFunction() {
            @Override
            public double evaluate(double argument) {
                double lnL = computeLogLikelihood(argument, taxonStates, nodeStates);

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

        double results[] = {minimum.minx / branchRates.getBranchRate(tree, node), -minimum.fminx};

        return results;
    }

    private AncestralStateBeagleTreeLikelihood asrLikelihood = null;
    private BranchRateModel branchRates = null;
    private PatternList patternList = null;
    private SubstitutionModel substitutionModel = null;
    private final Set<String> leafSet;
    private final Tree tree;
    private final DistanceType type;
    private final DataType dataType;
}
