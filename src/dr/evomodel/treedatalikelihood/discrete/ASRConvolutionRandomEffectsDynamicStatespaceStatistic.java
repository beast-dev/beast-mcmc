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

        double[] freqs = new double[64];
        for (int i = 0; i < 64; i++) {freqs[i] = 1.0/64.0;}
        this.tripletFreqs = new FrequencyModel(tripletDataType,freqs);

        this.tripletRates = new Parameter.Default("expandedRates",new double[4032]);
        this.tripletRatesRefx = new Parameter.Default("expandedRatesRefx",new double[4032]);
        this.tripletSubstitutionModel = new ComplexSubstitutionModel("expandedRateMatrix", tripletDataType, tripletFreqs, this.tripletRates);
        this.tripletSubstitutionModelRefx = new ComplexSubstitutionModel("expandedRateMatrixRefx", tripletDataType, tripletFreqs, this.tripletRates);

        // We normalize by hand to get triplets normalized per site
        tripletSubstitutionModel.setNormalization(false);
        tripletSubstitutionModelRefx.setNormalization(false);
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
        double[] freqsTriplet = new double[64];
        int idx = 0;
        double sum = 0.0;
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                for (int k = 0; k < 4; k++) {
                    freqsTriplet[idx] = freqsIID[i] * freqsIID[j] * freqsIID[k];
                    sum += freqsTriplet[idx];
                    idx++;
                }
            }
        }
        for (int i = 0; i < freqsTriplet.length; i++) {
            freqsTriplet[i] /= sum;
            tripletFreqs.setFrequency(i,freqsTriplet[i]);
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

    private int[] getTriplet(int i) {
        int[] from = new int[3];
        from[0] = Math.floorDiv(i, 16);
        from[1] = Math.floorDiv(i - (from[0] * 16), 4);
        from[2] = i - (from[0] * 16 + from[1] * 4);
        return from;
    }

    private void updateTransitionRateParameters(double[] matrixIID) {
        // TODO this is a really ugly function
        double[] freqsIID = subsModel.getFrequencyModel().getFrequencies();

        double sum1 = 0.0;
        double sum2 = 0.0;

        int idx = 0;
        // upper triangular
        for (int i = 0; i < 63; i++) {
            for (int j = i+1; j < 64; j++) {
                int[] from = getTriplet(i);
                int[] to = getTriplet(j);
                if (from[0] != to[0] && from[1] == to[1] && from[2] == to[2]) {
                    double bfProd = freqsIID[from[0]] * freqsIID[to[0]];
                    // infinitesimal matrix has stationary frequencies added but we don't want them
                    double rate = matrixIID[flatIndex(from[0], to[0])]/freqsIID[to[0]];
                    sum1 += rate * bfProd;
                    tripletRates.setParameterValue(i, rate);
                    rate *= Math.exp(randomEffects.getValue(i));
                    sum2 += rate * bfProd;
                    tripletRatesRefx.setParameterValue(i, rate);
                } else if (from[0] == to[0] && from[1] != to[1] && from[2] == to[2]) {
                    double bfProd = freqsIID[from[1]] * freqsIID[to[1]];
                    // infinitesimal matrix has stationary frequencies added but we don't want them
                    double rate = matrixIID[flatIndex(from[1], to[1])]/freqsIID[to[1]];
                    sum1 += rate * bfProd;
                    tripletRates.setParameterValue(i, rate);
                    rate *= Math.exp(randomEffects.getValue(i));
                    sum2 += rate * bfProd;
                    tripletRatesRefx.setParameterValue(i, rate);
                } else if (from[0] == to[0] && from[1] == to[1] && from[2] != to[2]) {
                    double bfProd = freqsIID[from[2]] * freqsIID[to[2]];
                    // infinitesimal matrix has stationary frequencies added but we don't want them
                    double rate = matrixIID[flatIndex(from[2], to[2])]/freqsIID[to[2]];
                    sum1 += rate * bfProd;
                    tripletRates.setParameterValue(i, rate);
                    rate *= Math.exp(randomEffects.getValue(i));
                    sum2 += rate * bfProd;
                    tripletRatesRefx.setParameterValue(i, rate);
                } else {
                    tripletRates.setParameterValue(i, 0.0);
                    tripletRatesRefx.setParameterValue(i, 0.0);
                }
                idx++;
            }
        }

        // lower triangular
        for (int i = 1; i < 64; i++) {
            for (int j = 0; j < i; j++) {
                int[] from = getTriplet(i);
                int[] to = getTriplet(j);
                if (from[0] != to[0] && from[1] == to[1] && from[2] == to[2]) {
                    double bfProd = freqsIID[from[0]] * freqsIID[to[0]];
                    // infinitesimal matrix has stationary frequencies added but we don't want them
                    double rate = matrixIID[flatIndex(from[0], to[0])]/freqsIID[to[0]];
                    sum1 += rate * bfProd;
                    tripletRates.setParameterValue(i, rate);
                    rate *= Math.exp(randomEffects.getValue(i));
                    sum2 += rate * bfProd;
                    tripletRatesRefx.setParameterValue(i, rate);
                } else if (from[0] == to[0] && from[1] != to[1] && from[2] == to[2]) {
                    double bfProd = freqsIID[from[1]] * freqsIID[to[1]];
                    // infinitesimal matrix has stationary frequencies added but we don't want them
                    double rate = matrixIID[flatIndex(from[1], to[1])]/freqsIID[to[1]];
                    sum1 += rate * bfProd;
                    tripletRates.setParameterValue(i, rate);
                    rate *= Math.exp(randomEffects.getValue(i));
                    sum2 += rate * bfProd;
                    tripletRatesRefx.setParameterValue(i, rate);
                } else if (from[0] == to[0] && from[1] == to[1] && from[2] != to[2]) {
                    double bfProd = freqsIID[from[2]] * freqsIID[to[2]];
                    // infinitesimal matrix has stationary frequencies added but we don't want them
                    double rate = matrixIID[flatIndex(from[2], to[2])]/freqsIID[to[2]];
                    sum1 += rate * bfProd;
                    tripletRates.setParameterValue(i, rate);
                    rate *= Math.exp(randomEffects.getValue(i));
                    sum2 += rate * bfProd;
                    tripletRatesRefx.setParameterValue(i, rate);
                } else {
                    tripletRates.setParameterValue(i, 0.0);
                    tripletRatesRefx.setParameterValue(i, 0.0);
                }
                idx++;
            }
        }

        // TODO is this the right normalization?
        for (int i = 0; i < 4032; i++) {
            tripletRates.setParameterValue(i, tripletRates.getParameterValue(i)/sum1);
            tripletRatesRefx.setParameterValue(i, tripletRatesRefx.getParameterValue(i)/sum2);
        }

    }

    private void updateSubstitutionRateMatrices() {
        // get triplet frequencies from site frequencies
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
        double[] tpm1 = new double[4096];
        double[] tpm2 = new double[4096];

        tripletSubstitutionModel.getTransitionProbabilities(distance1, tpm1);
        tripletSubstitutionModelRefx.getTransitionProbabilities(distance2, tpm2);

        // Do the matrix convolution
        double[] tpm= new double[4096];
        for (int i = 0; i < 64; i++) {
            for (int j = 0; j < 64; j++) {
                for (int k = 0; k < 64; k++) {
                    tpm[i * 64 + j] += tpm1[i * 64 + k] * tpm2[k * 64 + j];
                }
            }
        }

        double[] logTpm = new double[64 * 64];
        for (int i = 0; i < 64 * 64; i++) {
            logTpm[i] = Math.log(tpm[i]);
        }

        double lnL = 0.0;

        // We're assuming a low-mutation regime for this to work, so we assume ancestral and descendant "missing" nucleotide is the same
        // Then we average the transition probability over the missing state according to the per-site frequencies
        if ( ancestorStates[0] != descendantStates[0] ) {
            double pij = 0.0;
            for (int nuc = 0; nuc < 4; nuc++) {
                int ancestralTriplet = (nuc * 16) + (ancestorStates[0] * 4) + ancestorStates[1];
                int descendantTriplet = (nuc * 16) + (descendantStates[0] * 4) + descendantStates[1];
                pij += freqsIID[nuc] * tpm[ancestralTriplet * 64 + descendantTriplet];
            }
            lnL += Math.log(pij);
        }

        if ( ancestorStates[ancestorStates.length - 1] != descendantStates[ancestorStates.length - 1] ) {
            double pij = 0.0;
            int n = ancestorStates.length - 1;
            for (int nuc = 0; nuc < 4; nuc++) {
                int ancestralTriplet = (ancestorStates[n - 1] * 16) + (ancestorStates[n] * 4) + ancestorStates[nuc];
                int descendantTriplet = (ancestorStates[n - 1] * 16) + (ancestorStates[n] * 4) + descendantStates[nuc];
                pij += freqsIID[nuc] * tpm[ancestralTriplet * 64 + descendantTriplet];
            }
            lnL += Math.log(pij);
        }

        for (int s = 1; s < ancestorStates.length - 1; s++) {
            if ( ancestorStates[s] != descendantStates[s] ) {
                int ancestralTriplet = (ancestorStates[s-1] * 16) + (ancestorStates[s] * 4) + ancestorStates[s+1];
                int descendantTriplet = (descendantStates[s-1] * 16) + (descendantStates[s] * 4) + descendantStates[s+1];
                lnL += logTpm[ancestralTriplet * 64 + descendantTriplet];
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
    private final FrequencyModel tripletFreqs;
    private final Parameter randomEffects;
    private final Parameter tripletRates;
    private final Parameter tripletRatesRefx;
    private final ComplexSubstitutionModel tripletSubstitutionModel;
    private final ComplexSubstitutionModel tripletSubstitutionModelRefx;
    private final String[] sixtyFour = {"0","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42","43","44","45","46","47","48","49","50","51","52","53","54","55","56","57","58","59","60","61","62","63"};
    private final DataType tripletDataType = new GeneralDataType(sixtyFour);
}
