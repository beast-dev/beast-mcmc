package dr.evomodel.MSSD;

import dr.evolution.alignment.AscertainedSitePatterns;
import dr.evolution.alignment.PatternList;
import dr.evolution.datatype.MutationDeathType;
import dr.evolution.tree.NodeRef;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.sitemodel.SiteModel;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treelikelihood.LikelihoodCore;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.ParameterChangeType;
import dr.math.GammaFunction;

/**
 * Package: AbstractObservationProcess
 * Description:
 * <p/>
 * <p/>
 * Created by
 * Alexander V. Alekseyenko (alexander.alekseyenko@gmail.com)
 * Date: Feb 19, 2008
 * Time: 12:41:01 PM
 */
abstract public class AbstractObservationProcess extends AbstractModel {
    protected boolean[][] nodePatternInclusion;
    protected boolean[][] storedNodePatternInclusion;
    protected double[] cumLike;
    protected double[] nodePartials;
    protected double[] nodeLikelihoods;
    protected int nodeCount;
    protected int patternCount;
    protected int stateCount;
    protected TreeModel treeModel;
    protected PatternList patterns;
    protected double[] patternWeights;
    protected Parameter mu;
    protected Parameter lam;

    // update control variables
    protected boolean weightKnown;
    protected double logTreeWeight;
    protected double storedLogTreeWeight;
    private double gammaNorm;
    private double totalPatterns;
    protected MutationDeathType dataType;
    protected int deathState;
    protected SiteModel siteModel;
    private double logN;
    protected boolean nodePatternInclusionKnown = false;

    public AbstractObservationProcess(String Name, TreeModel treeModel, PatternList patterns, SiteModel siteModel,
                                      Parameter mu, Parameter lam) {
        super(Name);
        this.treeModel = treeModel;
        this.patterns = patterns;
        this.mu = mu;
        this.lam = lam;
        this.siteModel = siteModel;
        addModel(treeModel);
        addModel(siteModel);
        addParameter(mu);
        addParameter(lam);

        nodeCount = treeModel.getNodeCount();
        stateCount = patterns.getDataType().getStateCount();
        this.patterns = patterns;
        patternCount = patterns.getPatternCount();
        patternWeights = patterns.getPatternWeights();
        totalPatterns = 0;
        for (int i = 0; i < patternCount; ++i) {
            totalPatterns += patternWeights[i];
        }
        logN = Math.log(totalPatterns);

        gammaNorm = -GammaFunction.lnGamma(totalPatterns + 1);

        dataType = (MutationDeathType) patterns.getDataType();
        this.deathState = dataType.DEATHSTATE;
        setNodePatternInclusion();
        cumLike = new double[patternCount];
        nodeLikelihoods = new double[patternCount];
        weightKnown = false;
    }

    public Parameter getMuParameter() {
        return mu;
    }

    public Parameter getLamParameter() {
        return lam;
    }


    public final double nodePatternLikelihood(double[] freqs, BranchRateModel branchRateModel, LikelihoodCore likelihoodCore) {
        int i, j;
        double logL = gammaNorm;

        double birthRate = lam.getParameterValue(0);
        double prob;
        if (!nodePatternInclusionKnown)
            setNodePatternInclusion();
        if (nodePartials == null)
            nodePartials = new double[patternCount * stateCount];

        for (j = 0; j < patternCount; ++j) cumLike[j] = 0;

        for (i = 0; i < nodeCount; ++i) {
            // get partials for node i
            likelihoodCore.getPartials(i, nodePartials);
/*            System.err.println("Partials for node "+i);
            for(int m=0; m<patternCount;++m){
                System.err.print("[");
                for(int n=0; n<stateCount;++n){
                    System.err.print(" " +nodePartials[m*stateCount+n]);
                }
                System.err.println("]");
            }*/
            /*
                multiply the partials by equilibrium probs
                    this part could be optimized by first summing
                    and then multiplying by equilibrium probs
            */
            likelihoodCore.calculateLogLikelihoods(nodePartials, freqs, nodeLikelihoods);
            prob = Math.log(getNodeSurvivalProbability(i, branchRateModel));

            for (j = 0; j < patternCount; ++j) {
                if (nodePatternInclusion[i][j])
                    cumLike[j] += Math.exp(nodeLikelihoods[j] + prob);
            }
        }

        double ascertainmentCorrection = getAscertainmentCorrection(cumLike);
//        System.err.println("AscertainmentCorrection: "+ascertainmentCorrection);

        for (j = 0; j < patternCount; ++j) {
            logL += Math.log(cumLike[j] / ascertainmentCorrection) * patternWeights[j];
        }

        double deathRate = mu.getParameterValue(0);

        double logTreeWeight = getLogTreeWeight(branchRateModel);

/*        System.err.println("Patterns contribution "+logL);
        System.err.println("Patterns less gammanorm "+(logL-gammaNorm));
        System.err.println("TreeWeight "+getLogTreeWeight(branchRateModel));*/
        if (integrateGainRate) {
            logL -= gammaNorm + logN + Math.log(-logTreeWeight * deathRate / birthRate) * totalPatterns;
        } else {
            logL += logTreeWeight + Math.log(birthRate / deathRate) * totalPatterns;
        }
        return logL;
    }

    protected double getAscertainmentCorrection(double[] patternProbs) {
        // This function probably belongs better to the AscertainedSitePatterns
        double excludeProb = 0, includeProb = 0, returnProb = 1.0;
        if (this.patterns instanceof AscertainedSitePatterns) {
            int[] includeIndices = ((AscertainedSitePatterns) patterns).getIncludePatternIndices();
            int[] excludeIndices = ((AscertainedSitePatterns) patterns).getExcludePatternIndices();
            for (int i = 0; i < ((AscertainedSitePatterns) patterns).getIncludePatternCount(); i++) {
                int index = includeIndices[i];
                includeProb += patternProbs[index];
            }
            for (int j = 0; j < ((AscertainedSitePatterns) patterns).getExcludePatternCount(); j++) {
                int index = excludeIndices[j];
                excludeProb += patternProbs[index];
            }
            if (includeProb == 0.0) {
                returnProb -= excludeProb;
            } else if (excludeProb == 0.0) {
                returnProb = includeProb;
            } else {
                returnProb = includeProb - excludeProb;
            }
        }

        return returnProb;
    }

    final public double getLogTreeWeight(BranchRateModel branchRateModel) {
        if (!weightKnown) {
            logTreeWeight = calculateLogTreeWeight(branchRateModel);
            weightKnown = true;
        }

        return logTreeWeight;
    }

    abstract public double calculateLogTreeWeight(BranchRateModel branchRateModel);

    abstract void setNodePatternInclusion();

    final public double getAverageRate() {
        double avgRate = 0.0;
        double proportions[] = siteModel.getCategoryProportions();
        for (int i = 0; i < siteModel.getCategoryCount(); ++i) {
            avgRate += proportions[i] * siteModel.getRateForCategory(i);
        }
        return avgRate;
    }

    public double getNodeSurvivalProbability(int index, BranchRateModel branchRateModel) {
        NodeRef node = treeModel.getNode(index);
        NodeRef parent = treeModel.getParent(node);

        if (parent == null) return 1.0;

        final double deathRate = mu.getParameterValue(0) * getAverageRate();
        final double branchRate = branchRateModel.getBranchRate(treeModel, node);
        // Get the operational time of the branch
        final double branchTime = branchRate * treeModel.getBranchLength(node);
        return 1.0 - Math.exp(-deathRate * branchTime);
    }

    protected void handleModelChangedEvent(Model model, Object object, int index) {
        if (model == treeModel || model == siteModel) {
            weightKnown = false;
        }
        if (model == treeModel)
            nodePatternInclusionKnown = false;
    }

    protected final void handleParameterChangedEvent(Parameter parameter, int index, ParameterChangeType type) {
        if (parameter == mu || parameter == lam) {
            weightKnown = false;
        } else {
            System.err.println("AbstractObservationProcess: Got unexpected parameter changed event. (Parameter = " + parameter + ")");
        }
    }

    protected void storeState() {
        storedLogTreeWeight = logTreeWeight;
        storedNodePatternInclusion = nodePatternInclusion;
    }

    protected void restoreState() {
        logTreeWeight = storedLogTreeWeight;
        nodePatternInclusion = storedNodePatternInclusion;
    }

    protected void acceptState() {
    }

    public void setIntegrateGainRate(boolean integrateGainRate) {
        this.integrateGainRate = integrateGainRate;
    }

    private boolean integrateGainRate = false;
}
