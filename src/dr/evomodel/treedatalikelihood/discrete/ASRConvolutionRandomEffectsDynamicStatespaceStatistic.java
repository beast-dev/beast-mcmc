package dr.evomodel.treedatalikelihood.discrete;

import dr.evolution.datatype.DataType;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeUtils;
import dr.evolution.util.TaxonList;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evomodel.treelikelihood.AncestralStateBeagleTreeLikelihood;
import dr.inference.distribution.ParametricDistributionModel;
import dr.inference.model.Statistic;
import dr.math.MathUtils;
import dr.math.UnivariateFunction;
import dr.math.UnivariateMinimum;
import dr.xml.Reportable;

import java.util.Set;

/**
 * A statistic that returns the ML or MAP estimate of the times spent in two substitution models along a branch.
 *
 * @return The length of time before the descendant node at which the shift is estimated to occur.
 *
 * @author Andy Magee
 */
public class ASRConvolutionRandomEffectsDynamicStatespaceStatistic extends Statistic.Abstract implements Reportable {

    public ASRConvolutionRandomEffectsDynamicStatespaceStatistic(String name,
                                                                 AncestralStateBeagleTreeLikelihood asrLike,
                                                                 SubstitutionModel subsModelAncestor,
                                                                 SubstitutionModel subsModelDescendant,
                                                                 int[] pairFirstCharacters,
                                                                 int[] pairSecondCharacters,
                                                                 SubstitutionModel pairedSubsModelAncestor,
                                                                 SubstitutionModel pairedSubsModelDescendant,
                                                                 BranchRateModel branchRates,
                                                                 Statistic rateAncestor,
                                                                 Statistic rateDescendant,
                                                                 TaxonList mrcaTaxaDescendant,
                                                                 boolean bootstrap,
                                                                 ParametricDistributionModel prior
    ) throws TreeUtils.MissingTaxonException {
        this.name = name;
        this.bootstrap = bootstrap;
        this.prior = prior;
        this.asrLikelihood = asrLike;
        this.substitutionModelAncestor = subsModelAncestor;
        this.substitutionModelDescendant = subsModelDescendant;
        this.branchRates = branchRates;
        this.rateAncestor = rateAncestor;
        this.rateDescendant = rateDescendant;
        this.dataType = subsModelAncestor.getFrequencyModel().getDataType();
        if ( dataType != substitutionModelDescendant.getFrequencyModel().getDataType()) { throw new RuntimeException("Incompatible datatypes in substitution models for ASRSubstitutionModelConvolution.");}
        this.tree = asrLikelihood.getTreeModel();
        this.leafSetDescendant = (mrcaTaxaDescendant != null) ? TreeUtils.getLeavesForTaxa(tree, mrcaTaxaDescendant) : null;
        // Check validity of paired data model
        this.pairFirstCharacters = pairFirstCharacters;
        this.pairSecondCharacters = pairSecondCharacters;
        this.pairedSubstitutionModelAncestor = pairedSubsModelAncestor;
        this.pairedSubstitutionModelDescendant = pairedSubsModelDescendant;
        this.isPartitioned = pairedSubsModelAncestor != null;
        this.doubletsAreSafe = (!isPartitioned || !doubletsCanOverlap()) ? true : false;
        if (isPartitioned) {
            if ( bootstrap) { throw new RuntimeException("Cannot currently bootstrap context-dependent models."); }
            if ( (pairedSubsModelAncestor != null && pairedSubsModelAncestor == null) || pairedSubsModelAncestor == null && pairedSubsModelAncestor != null) { throw new RuntimeException("If specifying models for doublets must specify ancestral and descendant models."); }
            if ( pairedSubsModelAncestor.getFrequencyModel().getFrequencies().length != pairedSubsModelDescendant.getFrequencyModel().getFrequencies().length ) { throw new RuntimeException("Doublet models do not match in size."); }
            if ( pairedSubsModelAncestor.getFrequencyModel().getFrequencies().length != dataType.getStateCount() * dataType.getStateCount() ) { throw new RuntimeException("Doublet models are not sized for doublets."); }
        }
    }

    public int getDimension() {
        return 1;
    }

    public String getDimensionName(int i) {
        if ( name == null ) {
            return "timeBeforeMRCA";
        } else {
            return name;
        }

    }

    public String getStatisticName() {
        return NAME;
    }

    /**
     * @return the statistic
     */
    public double getStatisticValue(int dim) {
        NodeRef nodeDescendant = getNode(leafSetDescendant);
        NodeRef nodeAncestor = tree.getParent(nodeDescendant);

        if (!doubletsAreSafe && doubletsOverlapOnSequence(nodeAncestor)) {
            // We cannot technically apply this but to avoid mucking up a whole MCMC chain we'll just return NaN
            return Double.NaN;
        }

        double[] rates = getRates(nodeDescendant);

        double branchTime = tree.getNodeHeight(nodeAncestor) - tree.getNodeHeight(nodeDescendant);

        UnivariateMinimum optimized = optimizeTimes(nodeDescendant, nodeAncestor, rates, branchTime);

        return (1.0 - optimized.minx) * branchTime;
    }

    @Override
    public String getReport() {

        StringBuilder sb = new StringBuilder("hackyStatespaceStatistic Report\n\n");

        sb.append("Estimated time of shift before common ancestor: ").append(getStatisticValue(0)).append("\n");
        sb.append("Using rates: ").append(new dr.math.matrixAlgebra.Vector(getRates(getNode(leafSetDescendant)))).append("\n");
        sb.append("Using substitution models named: ").append(substitutionModelAncestor.getId()).append(", ").append(substitutionModelDescendant.getId()).append("\n");
        sb.append("Using taxon set: ").append(leafSetDescendant).append("\n");
        sb.append("Using prior? ").append((prior != null)).append("\n");
        if (prior != null) {
            sb.append("  Using prior of type: ").append((prior.getModelName())).append("\n");
            sb.append("  Using prior named: ").append((prior.getId())).append("\n");
        }
        sb.append("Using bootstrap? ").append((bootstrap)).append("\n");
        sb.append("\n\n");
        return sb.toString();
    }

    // Gets rates for ancestral and descendant (in that order) portions of the branch
    private double[] getRates(NodeRef node) {
        double[] rates = new double[2];

        if ( rateAncestor != null ) {
            rates[0] = rateAncestor.getStatisticValue(0);
        } else {
            rates[0] = branchRates.getBranchRate(tree,node);
        }

        if ( rateDescendant != null ) {
            rates[1] = rateDescendant.getStatisticValue(0);
        } else {
            rates[1] = branchRates.getBranchRate(tree,node);
        }

        return rates;
    }

    private void convolveMatrices(double[] tpm1, double[] tpm2, double[] tpm, int nStates) {
        for (int i = 0; i < nStates; i++) {
            for (int j = 0; j < nStates; j++) {
                for (int k = 0; k < nStates; k++) {
                    tpm[i * nStates + j] += tpm1[i * nStates + k] * tpm2[k * nStates + j];
                }
            }
        }
    }

    private boolean isPair(int first, int second) {
        if (!isPartitioned) {
            return false;
        }

        boolean pair = false;

        for (int i = 0; i < pairFirstCharacters.length; i++) {
            if ( first == pairFirstCharacters[i] && second == pairSecondCharacters[i] ) {
                pair = true;
                break;
            }
        }
        return pair;
    }

    private boolean doubletsCanOverlap() {
        boolean overlapPossible = false;
        for (int i = 0; i < pairFirstCharacters.length; i++) {
            for (int j= 0; j < pairSecondCharacters.length; j++) {
                if (pairFirstCharacters[i] == pairSecondCharacters[j]) {
                    overlapPossible = true;
                    break;
                }
            }
        }
        return overlapPossible;
    }

    private boolean doubletsOverlapOnSequence(NodeRef node) {
        int[] states = asrLikelihood.getStatesForNode(tree, node);
        boolean overlapping = false;
        boolean lastWasPair = false;
        for (int i = 0; i < states.length - 1; i++) {
            if ( isPair(states[i], states[i+1]) ) {
                if ( lastWasPair ) {
                    overlapping = true;
                    break;
                } else {
                    lastWasPair = true;
                }
            } else {
                lastWasPair = false;
            }
        }

        return overlapping;
    }

    private int getDoublet(int first, int second, int nStates) {
        return first * nStates + second;
    }

    private double computeLogLikelihood(double distance1, double distance2, double rateDescendant, int[] ancestorStates, int[] descendantStates) {
        int nStates = dataType.getStateCount();
        int nStatesSquared = nStates * nStates;

        double[] tpm1 = new double[nStatesSquared];
        double[] tpm2 = new double[nStatesSquared];
        substitutionModelAncestor.getTransitionProbabilities(distance1, tpm1);
        substitutionModelDescendant.getTransitionProbabilities(distance2, tpm2);

        // Do the matrix convolution
        double[] tpm= new double[nStates * nStates];
        convolveMatrices(tpm1, tpm2, tpm, nStates);

        double[] logTpm = new double[nStates * nStates];
        for (int i = 0; i < nStates * nStates; i++) {
            logTpm[i] = Math.log(tpm[i]);
        }

        double[] pairedTpm = new double[nStatesSquared * nStatesSquared];
        double[] logPairedTpm = new double[nStatesSquared * nStatesSquared];
        if ( isPartitioned ) {
            double[] pairedTpm1 = new double[nStatesSquared * nStatesSquared];
            double[] pairedTpm2 = new double[nStatesSquared * nStatesSquared];
            pairedSubstitutionModelAncestor.getTransitionProbabilities(distance1, pairedTpm1);
            pairedSubstitutionModelDescendant.getTransitionProbabilities(distance2, pairedTpm2);

            // Do the matrix convolution
            convolveMatrices(pairedTpm1, pairedTpm2, pairedTpm, nStates);

            for (int i = 0; i < nStatesSquared * nStatesSquared; i++) {
                logPairedTpm[i] = Math.log(pairedTpm1[i]);
            }
        }

        double lnL = 0.0;
        int nextToLast = ancestorStates.length - 1;
        for (int s = 0; s < ancestorStates.length; s++) {
            if (s < nextToLast && isPair(ancestorStates[s],ancestorStates[s+1])) {
                int from = getDoublet(ancestorStates[s], ancestorStates[s + 1], nStates);
                int to = getDoublet(descendantStates[s], descendantStates[s + 1], nStates);
                lnL += logPairedTpm[from * nStatesSquared + to];
            } else {
                lnL += logTpm[ancestorStates[s] * nStates + descendantStates[s]];
            }
        }

        if (prior != null) {
            // Prior is on (absolute) time of convolution before given MRCA
            lnL += prior.logPdf(distance2/rateDescendant);
        }

        return lnL;
    }

    private NodeRef getNode(Set<String> leafSet) {
        return (leafSet != null) ? TreeUtils.getCommonAncestorNode(tree, leafSet) : tree.getRoot();
    }

    private UnivariateMinimum optimizeTimes(NodeRef nodeDescendant, NodeRef nodeAncestor, double[] rates, double time) {
        int[] nodeStatesAncestor = asrLikelihood.getStatesForNode(tree, nodeAncestor);
        int[] nodeStatesDescendant = asrLikelihood.getStatesForNode(tree, nodeDescendant);

        int nStates = nodeStatesAncestor.length;

        if ( bootstrap ) {
            int[] tmp1 = new int[nStates];
            int[] tmp2 = new int[nStates];
            System.arraycopy(nodeStatesAncestor,0, tmp1,0, nStates);
            System.arraycopy(nodeStatesDescendant,0, tmp2,0, nStates);

            for (int i = 0; i < nStates; i++) {
                int idx = MathUtils.nextInt(nStates);
                nodeStatesAncestor[i] = tmp1[idx];
                nodeStatesDescendant[i] = tmp2[idx];
            }
        }

        UnivariateFunction f = new UnivariateFunction() {
            @Override
            public double evaluate(double argument) {
                double d1 = argument * rates[0] * time;
                double d2 = (1.0 - argument) * rates[1] * time;
                double lnL = computeLogLikelihood(d1, d2, rates[1], nodeStatesAncestor, nodeStatesDescendant);
//                System.err.println("Optimizing: time = " + time + ", rates[0] = " +
//                        rates[0] + ", rates[1] = " + rates[1] + ", proportion = " + argument +
//                        ", dist1 = " + d1 + ", dist2 = " + d2 + ", lnL = " + lnL);

                return -lnL;
            }

            @Override
            public double getLowerBound() {
                return 0;
            }

            @Override
            public double getUpperBound() {
                return 1.0;
            }
        };

        UnivariateMinimum minimum = new UnivariateMinimum();
        minimum.findMinimum(f);

        return minimum;
    }

    private BranchRateModel branchRates;
    private final Statistic rateAncestor;
    private final Statistic rateDescendant;
    private AncestralStateBeagleTreeLikelihood asrLikelihood;
    private SubstitutionModel substitutionModelAncestor;
    private SubstitutionModel substitutionModelDescendant;
    private SubstitutionModel pairedSubstitutionModelAncestor;
    private SubstitutionModel pairedSubstitutionModelDescendant;
    private final Set<String> leafSetDescendant;
    private final Tree tree;
    private final DataType dataType;
    private final boolean bootstrap;
    private final ParametricDistributionModel prior;
    private final String name;
    private final int[] pairFirstCharacters;
    private final int[] pairSecondCharacters;
    private final boolean isPartitioned;
    private final boolean doubletsAreSafe;
}
