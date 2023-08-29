package dr.evomodel.substmodel;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeUtils;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.coalescent.OldGMRFSkyrideLikelihood;
import dr.evomodel.siteratemodel.GammaSiteRateModel;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.tree.TreeStatistic;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.inference.markovjumps.MarkovJumpsCore;
import dr.inference.markovjumps.MarkovJumpsType;
import dr.inference.model.Parameter;
import dr.oldevomodel.sitemodel.GammaSiteModel;
import dr.oldevomodel.sitemodel.SiteModel;
import dr.xml.Reportable;

public class SubstitutionModelRandomEffectClassifier extends TreeStatistic implements Reportable {
    private final int dim;
    private final int nStates;
    private final int nSites;
    private final double threshold;

    private Tree tree;
    private GammaSiteRateModel siteModel;
    private BranchRateModel branchRates;

//    private ComplexSubstitutionModel proxy = null;
//    private Parameter relativeRateDummyParameter = null;
//    private MarkovJumpsSubstitutionModel markovJumps = null;
    
//    private double[] countMatrix;
    private int[] fromState;
    private int[] toState;

    public SubstitutionModelRandomEffectClassifier(String name,
                                                   Tree tree,
                                                   GammaSiteRateModel siteModel,
                                                   BranchRateModel branchRates,
                                                   int nSites,
                                                   double threshold) {
        super(name);
        this.tree = tree;
        this.siteModel = siteModel;
        this.branchRates = branchRates;

        if (!(siteModel.getSubstitutionModel() instanceof OldGLMSubstitutionModel)) {
            throw new RuntimeException("SubstitutionModelRandomEffectClassifier only works for GLM substitution models.");
        }

        this.nSites = nSites;
        this.nStates = siteModel.getSubstitutionModel().getFrequencyModel().getDataType().getStateCount();
        this.dim = nStates * (nStates - 1);
//        this.relativeRateDummyParameter = new Parameter.Default("",dim);
//        this.countMatrix = new double[nStates * nStates];
        this.threshold = threshold;

        this.fromState = new int[dim];
        this.toState = new int[dim];
        int idx = 0;
        int offset = dim/2;
        for (int i = 0; i < nStates - 1; i++) {
            for (int j = i + 1; j < nStates; j++) {
                fromState[idx] = i;
                toState[idx] = j;
                fromState[idx + offset] = j;
                toState[idx + offset] = i;
                idx++;
            }
        }
    }

    @Override
    public void setTree(Tree tree) {

    }

    @Override
    public Tree getTree() {
        return null;
    }

    @Override
    public int getDimension() {
        return dim;
    }

//    private void setupSubstitutionModels(int index, boolean includeRandomEffect) {
//        OldGLMSubstitutionModel glmSubs = (OldGLMSubstitutionModel)siteModel.getSubstitutionModel();
//        double[] relativeRates = new double[dim];
//        glmSubs.setupRelativeRates(relativeRates);
//
//        if (!includeRandomEffect) {
//            double[] copiedParameterValues = glmSubs.getGLM().getRandomEffect(0).getParameterValues();
////            double randomEffect = glmSubs.getGLM().getRandomEffect(0).getParameterValue(index);
//            relativeRates[index] /= Math.exp(copiedParameterValues[index]);
//        }
//
//        this.proxy = new ComplexSubstitutionModel("internalGlmProxyForSubstitutionModelRandomEffectClassifier",glmSubs.getDataType(),glmSubs.getFrequencyModel(),relativeRateDummyParameter);
//        this.markovJumps = new MarkovJumpsSubstitutionModel(proxy, MarkovJumpsType.COUNTS);
//    }

    private ComplexSubstitutionModel makeProxyModel(int index, boolean includeRandomEffect) {
        OldGLMSubstitutionModel glmSubs = (OldGLMSubstitutionModel)siteModel.getSubstitutionModel();
        double[] relativeRates = new double[dim];
        glmSubs.setupRelativeRates(relativeRates);

        if (!includeRandomEffect) {
            double[] copiedParameterValues = glmSubs.getGLM().getRandomEffect(0).getParameterValues();
//            double randomEffect = glmSubs.getGLM().getRandomEffect(0).getParameterValue(index);
            relativeRates[index] /= Math.exp(copiedParameterValues[index]);
        }

        Parameter relativeRateDummyParameter = new Parameter.Default(dim);
        for (int i = 0; i < dim; i++) {
            relativeRateDummyParameter.setParameterValueQuietly(i,relativeRates[i]);
        }

        return new ComplexSubstitutionModel("internalGlmProxyForSubstitutionModelRandomEffectClassifier",
                glmSubs.getDataType(),glmSubs.getFrequencyModel(),relativeRateDummyParameter);

    }

//    private MarkovJumpsSubstitutionModel makeMarkovJumpsSubstitutionModel(int index, boolean includeRandomEffect) {
//        return new MarkovJumpsSubstitutionModel(proxy, MarkovJumpsType.COUNTS);
//    }

//    private void setupCountMatrix(int index) {
//        for (int i = 0; i < dim; i++) {
//            countMatrix[i] = 0.0;
//        }
//        MarkovJumpsCore.fillRegistrationMatrix(countMatrix, fromState[index], toState[index], nStates, 1.0);
//    }
    private double[] makeCountMatrix(int index) {
        double[] countMatrix = new double[nStates * nStates];
        for (int i = 0; i < dim; i++) {
            countMatrix[i] = 0.0;
        }
        MarkovJumpsCore.fillRegistrationMatrix(countMatrix, fromState[index], toState[index], nStates, 1.0);
        return countMatrix;
    }

    private double getTreeLengthInSubstitutions(double relativeRate) {
        double length = 0.0;
        NodeRef root = tree.getRoot();
        for (int i = 0; i < tree.getNodeCount(); i++) {
            NodeRef node = tree.getNode(i);
            if ( node != root ) {
                double branchLength = tree.getNodeHeight(tree.getParent(node)) - tree.getNodeHeight(node);
                length += branchLength * branchRates.getBranchRate(tree, node);
            }
        }
        return length * relativeRate;
    }

    private double getDoubleResult(double countDiff) {
        if ( threshold <= 0.0 ) {
            return countDiff;
        } else {
            if ( Math.abs(countDiff) > threshold ) {
                return 1.0;
            } else {
                return 0.0;
            }
        }
    }

//    private double getCount(int index, boolean includeRandomEffect) {
//        setupSubstitutionModels(index, includeRandomEffect);
//        setupCountMatrix(index);
//        markovJumps.setRegistration(countMatrix);
//
//        double expectedCount = 0.0;
//        for (int i = 0; i < siteModel.getCategoryCount(); i++) {
//            double thisRateCategory = siteModel.getRateForCategory(i);
//            markovJumps.computeCondStatMarkovJumps(getTreeLengthInSubstitutions(thisRateCategory),countMatrix);
//            double countForThisCategory = countMatrix[fromState[index] * nStates + toState[index]];
//            expectedCount += countForThisCategory * siteModel.getProportionForCategory(i);
//        }
//
//        return expectedCount * nSites;
//    }
    
    private double getCountForCategory(int category, int index, boolean includeRandomEffect) {
        ComplexSubstitutionModel proxy = makeProxyModel(index, includeRandomEffect);
        MarkovJumpsSubstitutionModel markovJumps = new MarkovJumpsSubstitutionModel(proxy, MarkovJumpsType.COUNTS);

        double thisRateCategory = siteModel.getRateForCategory(category);
        double time = getTreeLengthInSubstitutions(thisRateCategory);

        double[] register = new double[nStates * nStates];
        double[] jointCounts = new double[nStates * nStates];
//        double[] conditionalCounts = new double[nStates * nStates];

        int from = fromState[index];
        int to = toState[index];
        MarkovJumpsCore.fillRegistrationMatrix(register, from, to, nStates, 1.0);
        markovJumps.setRegistration(register);

        markovJumps.computeJointStatMarkovJumps(time, jointCounts);

//        markovJumps.computeCondStatMarkovJumps(time, conditionalCounts);

        return jointCounts[from * nStates + to];
    }
    
    private double getCount(int index, boolean includeRandomEffect) {
        double expectedCount = 0.0;
        for (int i = 0; i < siteModel.getCategoryCount(); i++) {
            expectedCount += getCountForCategory(i, index, includeRandomEffect)  * siteModel.getProportionForCategory(i);
        }

        return expectedCount * nSites;
    }

    private double getCountDifferences(int index) {
        return getCount(index,true) - getCount(index, false);
    }

    @Override
    public double getStatisticValue(int dim) {
//        double[] countsWith = new double[this.dim];
//        double totCountsWith = 0.0;
//        double[] countsWithout = new double[this.dim];
//        double totCountsWithout = 0.0;
//        for (int i = 0; i < this.dim; i++) {
//            countsWith[i] = getCount(i,true);
//            countsWithout[i] = getCount(i,false);
//            totCountsWith += countsWith[i];
//            totCountsWithout += countsWithout[i];
//        }
//
//        System.err.println("Total counts without: " + totCountsWith + "\n");
//        System.err.println("Total counts with: " + totCountsWithout + "\n");
//
//        double meanTreeLength = 0.0;
//        for (int i = 0; i < siteModel.getCategoryCount(); i++) {
//            meanTreeLength += siteModel.getProportionForCategory(i) * getTreeLengthInSubstitutions(siteModel.getRateForCategory(i));
//        }

        return getDoubleResult(getCountDifferences(dim));
    }

}
