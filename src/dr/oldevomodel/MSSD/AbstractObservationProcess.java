/*
 * AbstractObservationProcess.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 *  BEAST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package dr.oldevomodel.MSSD;

import dr.evolution.alignment.AscertainedSitePatterns;
import dr.evolution.alignment.PatternList;
import dr.evolution.datatype.MutationDeathType;
import dr.evolution.tree.NodeRef;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.branchratemodel.DefaultBranchRateModel;
import dr.evomodel.tree.TreeChangedEvent;
import dr.oldevomodel.sitemodel.SiteRateModel;
import dr.evomodel.tree.TreeModel;
import dr.oldevomodel.treelikelihood.LikelihoodCore;
import dr.evomodel.treelikelihood.LikelihoodPartialsProvider;
import dr.oldevomodel.treelikelihood.ScaleFactorsHelper;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
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
    protected boolean[] nodePatternInclusion;
    protected boolean[] storedNodePatternInclusion;
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
    protected SiteRateModel siteModel;
    private double logN;
    protected boolean nodePatternInclusionKnown = false;
    BranchRateModel branchRateModel;

    public AbstractObservationProcess(String Name, TreeModel treeModel, PatternList patterns, SiteRateModel siteModel,
                                      BranchRateModel branchRateModel, Parameter mu, Parameter lam) {
        super(Name);
        this.treeModel = treeModel;
        this.patterns = patterns;
        this.mu = mu;
        this.lam = lam;
        this.siteModel = siteModel;
        if (branchRateModel != null) {
            this.branchRateModel = branchRateModel;
        } else {
            this.branchRateModel = new DefaultBranchRateModel();
        }
        addModel(treeModel);
        addModel(siteModel);
        addModel(this.branchRateModel);
        addVariable(mu);
        addVariable(lam);

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

//    public Parameter getMuParameter() {
//        return mu;
//    }
//
//    public Parameter getLamParameter() {
//        return lam;
//    }

    private double calculateSiteLogLikelihood(int site, double[] partials, double[] frequencies) {
        int v = site * stateCount;
        double sum = 0.0;
        for (int i = 0; i < stateCount; i++) {
            sum += frequencies[i] * partials[v + i];
        }
        return Math.log(sum);
    }


    private void calculateNodePatternLikelihood(int nodeIndex,
                                                double[] freqs,
                                                LikelihoodCore likelihoodCore,
                                                double averageRate,
                                                double[] cumLike) {
        // get partials for node nodeIndex
        likelihoodCore.getPartials(nodeIndex, nodePartials); // MAS
            /*
                multiply the partials by equilibrium probs
                    this part could be optimized by first summing
                    and then multiplying by equilibrium probs
            */
        double prob = Math.log(getNodeSurvivalProbability(nodeIndex, averageRate));

        for (int j = 0; j < patternCount; ++j) {
            if (nodePatternInclusion[nodeIndex * patternCount + j]) {
                cumLike[j] += Math.exp(calculateSiteLogLikelihood(j, nodePartials, freqs) + prob);
            }
        }
    }

    private double accumulateCorrectedLikelihoods(double[] cumLike, double ascertainmentCorrection,
                                                  double[] patterWeights) {
        double logL = 0;
        for (int j = 0; j < patternCount; ++j) {
            logL += Math.log(cumLike[j] / ascertainmentCorrection) * patternWeights[j];
        }
        return logL;
    }

    public final double nodePatternLikelihood(double[] freqs, LikelihoodPartialsProvider likelihoodCore,
                                              ScaleFactorsHelper scaleFactorsHelper) {
        int i, j;
        double logL = gammaNorm;

        double birthRate = lam.getParameterValue(0);
        double logProb;
        if (!nodePatternInclusionKnown)
            setNodePatternInclusion();
        if (nodePartials == null) {
            nodePartials = new double[patternCount * stateCount];
        }

        double averageRate = getAverageRate();

        for (j = 0; j < patternCount; ++j) cumLike[j] = 0;

        for (i = 0; i < nodeCount; ++i) {
            // get partials for node i
            likelihoodCore.getPartials(i, nodePartials);
            scaleFactorsHelper.rescalePartials(i, nodePartials);
            /*
                multiply the partials by equilibrium probs
                    this part could be optimized by first summing
                    and then multiplying by equilibrium probs
            */
//            likelihoodCore.calculateLogLikelihoods(nodePartials, freqs, nodeLikelihoods);   // MAS Removed
            logProb = Math.log(getNodeSurvivalProbability(i, averageRate));

            for (j = 0; j < patternCount; ++j) {
                if (nodePatternInclusion[i * patternCount + j]) {
//                    cumLike[j] += Math.exp(nodeLikelihoods[j] + logProb);  // MAS Replaced with line below
                    cumLike[j] += Math.exp(calculateSiteLogLikelihood(j, nodePartials, freqs)
                            + logProb);
                }
            }
        }

        double ascertainmentCorrection = getAscertainmentCorrection(cumLike);
//        System.err.println("AscertainmentCorrection: "+ascertainmentCorrection);

        for (j = 0; j < patternCount; ++j) {
            logL += Math.log(cumLike[j] / ascertainmentCorrection) * patternWeights[j];
        }

        double deathRate = mu.getParameterValue(0);

        double logTreeWeight = getLogTreeWeight();

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

    final public double getLogTreeWeight() {
        if (!weightKnown) {
            logTreeWeight = calculateLogTreeWeight();
            weightKnown = true;
        }

        return logTreeWeight;
    }

    abstract public double calculateLogTreeWeight();

    abstract void setNodePatternInclusion();

    final public double getAverageRate() {
        if (!averageRateKnown) {
            double avgRate = 0.0;
            double proportions[] = siteModel.getCategoryProportions();
            for (int i = 0; i < siteModel.getCategoryCount(); ++i) {
                avgRate += proportions[i] * siteModel.getRateForCategory(i);
            }
            averageRate = avgRate;
            averageRateKnown = true;
        }
        return averageRate;
    }

    public double getNodeSurvivalProbability(int index, double averageRate) {
        NodeRef node = treeModel.getNode(index);
        NodeRef parent = treeModel.getParent(node);

        if (parent == null) return 1.0;

        final double deathRate = mu.getParameterValue(0) * averageRate; //getAverageRate();
        final double branchRate = branchRateModel.getBranchRate(treeModel, node);
        // Get the operational time of the branch
        final double branchTime = branchRate * treeModel.getBranchLength(node);
        return 1.0 - Math.exp(-deathRate * branchTime);
    }

    protected void handleModelChangedEvent(Model model, Object object, int index) {
        if (model == siteModel) {
            averageRateKnown = false;
        }
        if (model == treeModel || model == siteModel || model == branchRateModel) {
            weightKnown = false;
        }
        if (model == treeModel) {
            if (object instanceof TreeChangedEvent) {
                if (((TreeChangedEvent) object).isTreeChanged()) {
                    nodePatternInclusionKnown = false;
                }
            }
        }
    }

    protected final void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        if (variable == mu || variable == lam) {
            weightKnown = false;
        } else {
            System.err.println("AbstractObservationProcess: Got unexpected parameter changed event. (Parameter = " + variable + ")");
        }
    }

    protected void storeState() {
//        storedAverageRate = averageRate;
        storedLogTreeWeight = logTreeWeight;
        System.arraycopy(nodePatternInclusion, 0, storedNodePatternInclusion, 0, storedNodePatternInclusion.length);
    }

    protected void restoreState() {
//        averageRate = storedAverageRate;
        averageRateKnown = false;
        logTreeWeight = storedLogTreeWeight;
        boolean[] tmp = storedNodePatternInclusion;
        storedNodePatternInclusion = nodePatternInclusion;
        nodePatternInclusion = tmp;
    }

    protected void acceptState() {
    }

    public void setIntegrateGainRate(boolean integrateGainRate) {
        this.integrateGainRate = integrateGainRate;
    }

    private boolean integrateGainRate = false;

    private double storedAverageRate;
    private double averageRate;
    private boolean averageRateKnown = false;
}
