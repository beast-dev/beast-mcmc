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
                                                    SubstitutionModel subsModelAncestor,
                                                    SubstitutionModel subsModelDescendant,
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
        this.subsModelAncestor = subsModelAncestor;
        this.subsModelDescendant = subsModelDescendant;
        this.randomEffects = randomEffects;
        this.branchRates = branchRates;
        this.rateAncestor = rateAncestor;
        this.rateDescendant = rateDescendant;
        this.tree = asrLikelihood.getTreeModel();
        this.leafSetDescendant = (mrcaTaxaDescendant != null) ? TreeUtils.getLeavesForTaxa(tree, mrcaTaxaDescendant) : null;

        if (!(subsModelAncestor.getFrequencyModel().getDataType() instanceof Nucleotides)) {throw new RuntimeException("This statistic only works for DNA data");}
        if (!(subsModelDescendant.getFrequencyModel().getDataType() instanceof Nucleotides)) {throw new RuntimeException("This statistic only works for DNA data");}

        double[] freqs = new double[16];
        for (int i = 0; i < 16; i++) {freqs[i] = 1.0/16.0;}
        this.doubletFreqs = new FrequencyModel(doubletDataType,freqs);

        this.doubletRates = new Parameter.Default("expandedRates",new double[240]);
        this.doubletRatesRefx = new Parameter.Default("expandedRatesRefx",new double[240]);
        this.doubletSubstitutionModel = new ComplexSubstitutionModel("expandedRateMatrix", doubletDataType, doubletFreqs, doubletRates);
        this.doubletSubstitutionModelRefx = new ComplexSubstitutionModel("expandedRateMatrixRefx", doubletDataType, doubletFreqs, doubletRatesRefx);

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

        updateSubstitutionRateMatrices();
        double[] tmp = new double[256];
        doubletSubstitutionModel.getInfinitesimalMatrix(tmp);
//        System.err.println(new dr.math.matrixAlgebra.Vector(tmp));
//        System.err.println(new dr.math.matrixAlgebra.Matrix(tmp,16,16));
        doubletSubstitutionModelRefx.getInfinitesimalMatrix(tmp);
//        System.err.println(new dr.math.matrixAlgebra.Vector(tmp));
//        System.err.println(new dr.math.matrixAlgebra.Matrix(tmp,16,16));

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

    private int[] nucsFromDoublet(int i) {
        int[] nucs = new int[2];
        nucs[0] = Math.floorDiv(i, 4);
        nucs[1] = i - (nucs[0] * 4);
        return nucs;
    }

    private int doubletFromNucs(int nuc1, int nuc2) {
        return 4 * nuc1 + nuc2;
    }

    private void updateTransitionRateParameters() {
        // TODO this is a really ugly function
        double[] matrixIID = new double[16];
        subsModelAncestor.getInfinitesimalMatrix(matrixIID);

//        System.err.println(new dr.math.matrixAlgebra.Vector(matrixIID));
//        System.err.println("\n\n\n");
//        System.err.println(new dr.math.matrixAlgebra.Matrix(matrixIID,4,4));

        double[] freqsIID = subsModelAncestor.getFrequencyModel().getFrequencies();
        double[] freqs = doubletFreqs.getFrequencies();

        double sum1 = 0.0;
        double sum2 = 0.0;

        int idx = 0;
        // upper triangular
        for (int i = 0; i < 15; i++) {
            for (int j = i+1; j < 16; j++) {
                int[] from = nucsFromDoublet(i);
                int[] to = nucsFromDoublet(j);
                if (from[0] != to[0] && from[1] == to[1]) {
                    double bfProd = freqs[from[0]] * freqs[to[0]];
                    // infinitesimal matrix has stationary frequencies added but we don't want them
                    double rate = matrixIID[flatIndex(from[0], to[0])]/freqsIID[to[0]];
//                    System.err.println("  (0) rate = " + rate + "; idx = " + idx + "; refx = " + randomEffects.getParameterValue(idx));
                    sum1 += rate * bfProd;
                    doubletRates.setParameterValue(idx, rate);
                    rate *= Math.exp(randomEffects.getParameterValue(idx));
                    sum2 += rate * bfProd;
                    doubletRatesRefx.setParameterValue(idx, rate);
//                    System.err.println("rate * refx = " + rate);
                } else if (from[0] == to[0] && from[1] != to[1]) {
                    double bfProd = freqs[from[0]] * freqs[to[0]];
                    // infinitesimal matrix has stationary frequencies added but we don't want them
                    double rate = matrixIID[flatIndex(from[1], to[1])]/freqsIID[to[1]];
//                    System.err.println("  (1) rate = " + rate + "; idx = " + idx + "; refx = " + randomEffects.getParameterValue(idx));
                    sum1 += rate * bfProd;
                    doubletRates.setParameterValue(idx, rate);
                    rate *= Math.exp(randomEffects.getParameterValue(idx));
                    sum2 += rate * bfProd;
                    doubletRatesRefx.setParameterValue(idx, rate);
//                    System.err.println("rate * refx = " + rate);
                } else {
                    doubletRates.setParameterValue(idx, 0.0);
                    doubletRatesRefx.setParameterValue(idx, 0.0);
                }
                idx++;
            }
        }
//        System.err.println("\n>>>LOWER<<<\n");
        // lower triangular
        for (int j = 0; j < 15; j++) {
            for (int i = j+1; i < 16; i++) {
                int[] from = nucsFromDoublet(i);
                int[] to = nucsFromDoublet(j);
//                System.err.println("idx = " + idx +"; i = " + i + "; j = " + j + "; from = (" + from[0] + "," + from[1] + "); to = (" + to[0] + "," + to[1] + ")");
                if (from[0] != to[0] && from[1] == to[1]) {
                    double bfProd = freqs[from[0]] * freqs[to[0]];
                    // infinitesimal matrix has stationary frequencies added but we don't want them
                    double rate = matrixIID[flatIndex(from[0], to[0])]/freqsIID[to[0]];
//                    System.err.println("  (0) rate = " + rate + "; idx = " + idx + "; refx = " + randomEffects.getParameterValue(idx));
                    sum1 += rate * bfProd;
                    doubletRates.setParameterValue(idx, rate);
                    rate *= Math.exp(randomEffects.getParameterValue(idx));
                    sum2 += rate * bfProd;
                    doubletRatesRefx.setParameterValue(idx, rate);
//                    System.err.println("rate * refx = " + rate);
                } else if (from[0] == to[0] && from[1] != to[1]) {
                    double bfProd = freqs[from[0]] * freqs[to[0]];
                    // infinitesimal matrix has stationary frequencies added but we don't want them
                    double rate = matrixIID[flatIndex(from[1], to[1])]/freqsIID[to[1]];
//                    System.err.println("  (1) rate = " + rate + "; idx = " + idx + "; refx = " + randomEffects.getParameterValue(idx));
                    sum1 += rate * bfProd;
                    doubletRates.setParameterValue(idx, rate);
                    rate *= Math.exp(randomEffects.getParameterValue(idx));
                    sum2 += rate * bfProd;
                    doubletRatesRefx.setParameterValue(idx, rate);
//                    System.err.println("rate * refx = " + rate);
                } else {
                    doubletRates.setParameterValue(idx, 0.0);
                    doubletRatesRefx.setParameterValue(idx, 0.0);
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
        updateFrequencyModel(subsModelAncestor.getFrequencyModel().getFrequencies());
        updateTransitionRateParameters();
    }

    private void doConvolution(double[] tpm1, double[] tpm2, double[] tpm, int nStates) {
        for (int i = 0; i < nStates; i++) {
            for (int j = 0; j < nStates; j++) {
                for (int k = 0; k < nStates; k++) {
                    tpm[i * nStates + j] += tpm1[i * nStates + k] * tpm2[k * nStates + j];
                }
            }
        }
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

        double[] freqsIID = subsModelAncestor.getFrequencyModel().getFrequencies();

        double[] tpmDoublet1 = new double[256];
        double[] tpmDoublet2 = new double[256];
        doubletSubstitutionModel.getTransitionProbabilities(distance1, tpmDoublet1);
        doubletSubstitutionModelRefx.getTransitionProbabilities(distance2, tpmDoublet2);

        double[] tpmPerSite1 = new double[16];
        double[] tpmPerSite2 = new double[16];
        subsModelAncestor.getTransitionProbabilities(distance1, tpmPerSite1);
        subsModelDescendant.getTransitionProbabilities(distance2, tpmPerSite2);

        // Do the matrix convolution
        double[] tpmDoublet = new double[256];
        double[] tpmPerSite = new double[16];

        doConvolution(tpmDoublet1, tpmDoublet2, tpmDoublet, 16);
        doConvolution(tpmPerSite1, tpmPerSite2, tpmPerSite, 4);


        double[] logTpmPerSite = new double[16];
        for (int i = 0; i < 16; i++) {
            logTpmPerSite[i] = Math.log(tpmPerSite[i]);
        }

        double[] logTpmDoublet = new double[256];
        for (int i = 0; i < 256; i++) {
            logTpmDoublet[i] = Math.log(tpmDoublet[i]);
        }

        double lnL = 0.0;

        // TODO always check that we don't have 2 in a row mutations first!!!!

        double lastLnL = logTpmPerSite[4 * ancestorStates[0] + descendantStates[0]];
        lnL += lastLnL;
        for (int s = 1; s < ancestorStates.length; s++) {
            if ( ancestorStates[s] != descendantStates[s] ) {
                // Assign GA pairs to the dinucleotide model
                if (ancestorStates[s] == 2 && ancestorStates[s + 1] == 0) {
                    int ancestralDoublet = doubletFromNucs(ancestorStates[s], ancestorStates[s + 1]);
                    int descendantDoublet = doubletFromNucs(descendantStates[s], descendantStates[s + 1]);
                    lastLnL = logTpmDoublet[16 * ancestralDoublet + descendantDoublet];
                    // Skip next site, we've already computed its likelihood
                    s++;
                // Assign TT pairs to the dinucleotide model
                } else if (ancestorStates[s] == 3 && ancestorStates[s - 1] == 3) {
                    // Remove last site and add it to the doublet model
                    lnL -= lastLnL;
                    int ancestralDoublet = doubletFromNucs(ancestorStates[s - 1], ancestorStates[s]);
                    int descendantDoublet = doubletFromNucs(descendantStates[s - 1], descendantStates[s]);
                    lastLnL = logTpmDoublet[16 * ancestralDoublet + descendantDoublet];
                // Everything else gets the single site model
                } else {
                    lastLnL = logTpmPerSite[4 * ancestorStates[s] + descendantStates[s]];
                }
            } else {
                lastLnL = logTpmPerSite[4 * ancestorStates[s] + descendantStates[s]];
            }
            lnL += lastLnL;
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
    private SubstitutionModel subsModelAncestor;
    private SubstitutionModel subsModelDescendant;
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
