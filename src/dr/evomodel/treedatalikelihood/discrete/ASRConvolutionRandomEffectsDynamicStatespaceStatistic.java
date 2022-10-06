package dr.evomodel.treedatalikelihood.discrete;

        import dr.evolution.datatype.Codons;
        import dr.evolution.datatype.DataType;
        import dr.evolution.datatype.GeneralDataType;
        import dr.evolution.datatype.Nucleotides;
        import dr.evolution.tree.NodeRef;
        import dr.evolution.tree.Tree;
        import dr.evolution.tree.TreeUtils;
        import dr.evolution.util.TaxonList;
        import dr.evomodel.branchratemodel.BranchRateModel;
        import dr.evomodel.substmodel.*;
        import dr.evomodel.treelikelihood.AncestralStateBeagleTreeLikelihood;
        import dr.inference.distribution.ParametricDistributionModel;
        import dr.inference.glm.GeneralizedLinearModel;
        import dr.inference.hmc.CompoundGradient;
        import dr.inference.model.Parameter;
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
                                                    SubstitutionModel subsModel,
                                                    Parameter randomEffects,
                                                    BranchRateModel branchRates,
                                                    Statistic rateAncestor,
                                                    Statistic rateDescendant,
                                                    TaxonList mrcaTaxaDescendant,
                                                    ParametricDistributionModel prior
    ) throws TreeUtils.MissingTaxonException {
        this.name = name;
        this.prior = prior;
        this.asrLikelihood = asrLike;
        this.subsModel = subsModel;
        this.randomEffects = randomEffects;
        this.branchRates = branchRates;
        this.rateAncestor = rateAncestor;
        this.rateDescendant = rateDescendant;
        this.tree = asrLikelihood.getTreeModel();
        this.leafSetDescendant = (mrcaTaxaDescendant != null) ? TreeUtils.getLeavesForTaxa(tree, mrcaTaxaDescendant) : null;

        if (!(subsModel.getFrequencyModel().getDataType() instanceof Nucleotides)) {throw new RuntimeException("This statistic only works for DNA data");}

        double[] freqs = new double[16];
        for (int i = 0; i < 16; i++) {freqs[i] = 1.0/16.0;}
        this.doubletFreqs = new FrequencyModel(doubletDataType,freqs);

        this.doubletRates = new Parameter.Default("expandedRates",new double[240]);
        this.doubletRatesRefx = new Parameter.Default("expandedRatesRefx",new double[240]);
        this.doubletSubstitutionModel = new ComplexSubstitutionModel("expandedRateMatrix", doubletDataType, doubletFreqs, this.doubletRates);
        this.doubletSubstitutionModelRefx = new ComplexSubstitutionModel("expandedRateMatrixRefx", doubletDataType, doubletFreqs, this.doubletRates);

        // We normalize by hand to get doublets normalized per site
        doubletSubstitutionModel.setNormalization(false);
        doubletSubstitutionModelRefx.setNormalization(false);
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

        double[] rates = getRates(nodeDescendant);

        double branchTime = tree.getNodeHeight(nodeAncestor) - tree.getNodeHeight(nodeDescendant);

        UnivariateMinimum optimized = optimizeTimes(nodeDescendant, nodeAncestor, rates, branchTime);

        return (1.0 - optimized.minx) * branchTime;
    }

    @Override
    public String getReport() {

        StringBuilder sb = new StringBuilder("ASRSubstitutionModelConvolutionStatistic Report\n\n");

        sb.append("Estimated time of shift before common ancestor: ").append(getStatisticValue(0)).append("\n");
        sb.append("Using rates: ").append(new dr.math.matrixAlgebra.Vector(getRates(getNode(leafSetDescendant)))).append("\n");
        sb.append("Using taxon set: ").append(leafSetDescendant).append("\n");
        sb.append("Using prior? ").append((prior != null)).append("\n");
        if (prior != null) {
            sb.append("  Using prior of type: ").append((prior.getModelName())).append("\n");
            sb.append("  Using prior named: ").append((prior.getId())).append("\n");
        }
        sb.append("\n\n");
        return sb.toString();
    }

    private void updateFrequencyModel(double[] freqsIID) {
        double[] freqsDoublet = new double[16];
        int idx = 0;
        double sum = 0.0;
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                freqsDoublet[idx] = freqsIID[i] * freqsIID[j];
                sum += freqsDoublet[idx];
                idx++;
            }
        }
        for (int i = 0; i < freqsDoublet.length; i++) {
            freqsDoublet[i] /= sum;
            doubletFreqs.setFrequency(i,freqsDoublet[i]);
        }
    }

//    private int flatIndex(int i, int j) {
//        int k = 0;
//        if (i > j) {
//            k += 6;
//            int tmp = i;
//            i = j;
//            j = tmp;
//        }
//        k += (4*(4-1)/2) - (4-i)*((4-i)-1)/2 + j - i - 1;
//        return k;
//    }

    private int flatIndex(int i, int j) {
        return 4 * i + j;
    }

    private int[] getDoublet(int i) {
        int[] from = new int[2];
        from[0] = Math.floorDiv(i, 4);
        from[1] = i - (from[0] * 4);
        return from;
    }

    private void updateTransitionRateParameters(double[] matrixIID) {
        // TODO this is a really ugly function
        double[] freqsIID = subsModel.getFrequencyModel().getFrequencies();
        double[] freqs = doubletFreqs.getFrequencies();

        double sum1 = 0.0;
        double sum2 = 0.0;

        int idx = 0;
        // upper triangular
        for (int i = 0; i < 15; i++) {
            for (int j = i+1; j < 16; j++) {
                int[] from = getDoublet(i);
                int[] to = getDoublet(j);
                if (from[0] != to[0] && from[1] == to[1]) {
                    double bfProd = freqs[from[0]] * freqsIID[to[0]];
                    // infinitesimal matrix has stationary frequencies added but we don't want them
                    double rate = matrixIID[flatIndex(from[0], to[0])]/freqsIID[to[0]];
                    sum1 += rate * bfProd;
                    doubletRates.setParameterValue(i, rate);
                    rate *= Math.exp(randomEffects.getValue(i));
                    sum2 += rate * bfProd;
                    doubletRatesRefx.setParameterValue(i, rate);
                } else if (from[0] == to[0] && from[1] != to[1]) {
                    double bfProd = freqsIID[from[1]] * freqsIID[to[1]];
                    // infinitesimal matrix has stationary frequencies added but we don't want them
                    double rate = matrixIID[flatIndex(from[1], to[1])]/freqsIID[to[1]];
                    sum1 += rate * bfProd;
                    doubletRates.setParameterValue(i, rate);
                    rate *= Math.exp(randomEffects.getValue(i));
                    sum2 += rate * bfProd;
                    doubletRatesRefx.setParameterValue(i, rate);
                } else {
                    doubletRates.setParameterValue(i, 0.0);
                    doubletRatesRefx.setParameterValue(i, 0.0);
                }
                idx++;
            }
        }

        // lower triangular
        for (int i = 1; i < 16; i++) {
            for (int j = 0; j < i; j++) {
                int[] from = getDoublet(i);
                int[] to = getDoublet(j);
                if (from[0] != to[0] && from[1] == to[1]) {
                    double bfProd = freqs[from[0]] * freqsIID[to[0]];
                    // infinitesimal matrix has stationary frequencies added but we don't want them
                    double rate = matrixIID[flatIndex(from[0], to[0])]/freqsIID[to[0]];
                    sum1 += rate * bfProd;
                    doubletRates.setParameterValue(i, rate);
                    rate *= Math.exp(randomEffects.getValue(i));
                    sum2 += rate * bfProd;
                    doubletRatesRefx.setParameterValue(i, rate);
                } else if (from[0] == to[0] && from[1] != to[1]) {
                    double bfProd = freqsIID[from[1]] * freqsIID[to[1]];
                    // infinitesimal matrix has stationary frequencies added but we don't want them
                    double rate = matrixIID[flatIndex(from[1], to[1])]/freqsIID[to[1]];
                    sum1 += rate * bfProd;
                    doubletRates.setParameterValue(i, rate);
                    rate *= Math.exp(randomEffects.getValue(i));
                    sum2 += rate * bfProd;
                    doubletRatesRefx.setParameterValue(i, rate);
                } else {
                    doubletRates.setParameterValue(i, 0.0);
                    doubletRatesRefx.setParameterValue(i, 0.0);
                }
                idx++;
            }
        }

        // TODO is this the right normalization?
        for (int i = 0; i < 240; i++) {
            doubletRates.setParameterValue(i, doubletRates.getParameterValue(i)/sum1);
            doubletRatesRefx.setParameterValue(i, doubletRatesRefx.getParameterValue(i)/sum2);
        }

    }

    private void updateSubstitutionRateMatrices() {
        // get doublet frequencies from site frequencies
        updateFrequencyModel(subsModel.getFrequencyModel().getFrequencies());

        double[] rateMatrixIIDFlat = new double[16];
        subsModel.getInfinitesimalMatrix(rateMatrixIIDFlat);
        updateTransitionRateParameters(rateMatrixIIDFlat);
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

    private double computeLogLikelihood(double distance1, double distance2, double rateDescendant, int[] ancestorStates, int[] descendantStates) {

        double[] freqsIID = subsModel.getFrequencyModel().getFrequencies();

        updateSubstitutionRateMatrices();
        double[] tpm1 = new double[256];
        double[] tpm2 = new double[256];

        doubletSubstitutionModel.getTransitionProbabilities(distance1, tpm1);
        doubletSubstitutionModelRefx.getTransitionProbabilities(distance2, tpm2);

        // Do the matrix convolution
        double[] tpm= new double[256];
        for (int i = 0; i < 16; i++) {
            for (int j = 0; j < 16; j++) {
                for (int k = 0; k < 16; k++) {
                    tpm[i * 16 + j] += tpm1[i * 16 + k] * tpm2[k * 16 + j];
                }
            }
        }

        double[] logTpm = new double[16 * 16];
        for (int i = 0; i < 16 * 16; i++) {
            logTpm[i] = Math.log(tpm[i]);
        }

        double lnL = 0.0;

        // We're assuming a low-mutation regime for this to work, so we assume ancestral and descendant "missing" nucleotide is the same
        // Then we average the transition probability over the missing state according to the per-site frequencies
        if ( ancestorStates[0] != descendantStates[0] ) {
            double pij = 0.0;
            for (int nuc = 0; nuc < 4; nuc++) {
                int ancestralDoublet = (nuc * 4) + ancestorStates[0];
                int descendantDoublet = (nuc * 4) + descendantStates[0];
                pij += freqsIID[nuc] * tpm[ancestralDoublet * 16 + descendantDoublet];
            }
            lnL += Math.log(pij);
        }

        for (int s = 1; s < ancestorStates.length; s++) {
            if ( ancestorStates[s] != descendantStates[s] ) {
                int ancestralDoublet = (ancestorStates[s-1] * 4) + ancestorStates[s];
                int descendantDoublet = (descendantStates[s-1] * 4) + descendantStates[s];
                lnL += logTpm[ancestralDoublet * 16 + descendantDoublet];
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
    private SubstitutionModel subsModel;
    private final Set<String> leafSetDescendant;
    private final Tree tree;
    private final ParametricDistributionModel prior;
    private final String name;
    private final FrequencyModel doubletFreqs;
    private final Parameter randomEffects;
    private final Parameter doubletRates;
    private final Parameter doubletRatesRefx;
    private final ComplexSubstitutionModel doubletSubstitutionModel;
    private final ComplexSubstitutionModel doubletSubstitutionModelRefx;
    private final String[] sixteen = {"0","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15"};
    private final DataType doubletDataType = new GeneralDataType(sixteen);
}
