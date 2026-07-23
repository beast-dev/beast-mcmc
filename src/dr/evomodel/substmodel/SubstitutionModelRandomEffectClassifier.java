/*
 * SubstitutionModelRandomEffectClassifier.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
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
 *
 */

package dr.evomodel.substmodel;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.branchmodel.BranchModel;
import dr.evomodel.branchmodel.EpochBranchModel;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.siteratemodel.GammaSiteRateModel;
import dr.evomodel.tree.TreeStatistic;
import dr.inference.markovjumps.MarkovJumpsCore;
import dr.inference.markovjumps.MarkovJumpsType;
import dr.inference.model.Parameter;
import dr.xml.Reportable;

import java.util.List;

public class SubstitutionModelRandomEffectClassifier extends TreeStatistic implements Reportable {
    private final int dim;
    private final int nStates;
    private final int nSites;
    private final double threshold;

    private Tree tree;
    private GlmSubstitutionModel glmSubstitutionModel;
    private EpochBranchModel epochBranchModel;
    private GammaSiteRateModel siteModel;
    private BranchRateModel branchRates;

    private Parameter proxyRates;
    private ComplexSubstitutionModel proxy;
    private MarkovJumpsSubstitutionModel markovJumps;

    private final boolean usingRateVariation;
    private final boolean usingEpochs;
    private boolean[] epochUsesTargetModel;
    private boolean nullIsZero;
//    private double[] countMatrix;
    private int[] fromState;
    private int[] toState;

    public SubstitutionModelRandomEffectClassifier(String name,
                                                   Tree tree,
                                                   GlmSubstitutionModel glmSubstitutionModel,
                                                   EpochBranchModel epochBranchModel,
                                                   BranchRateModel branchRates,
                                                   GammaSiteRateModel siteModel,
                                                   int nSites,
                                                   double threshold,
                                                   boolean nullIsZero) {
        super(name);
        this.tree = tree;
        this.glmSubstitutionModel = glmSubstitutionModel;
        this.epochBranchModel = epochBranchModel;
        this.siteModel = siteModel;
        this.branchRates = branchRates;
        this.nullIsZero = nullIsZero;

        usingRateVariation = siteModel != null ? true : false;
        usingEpochs = epochBranchModel != null ? true : false;

        if (!(glmSubstitutionModel instanceof GlmSubstitutionModel)) {
            throw new RuntimeException("SubstitutionModelRandomEffectClassifier only works for GLM substitution models.");
        }

        if (usingEpochs) {
            List<SubstitutionModel> substitutionModels = epochBranchModel.getSubstitutionModels();
            this.epochUsesTargetModel = new boolean[substitutionModels.size()];
            int matches = 0;
            for (int i = 0; i < substitutionModels.size(); i++) {
                if (substitutionModels.get(i) == glmSubstitutionModel) {
                    this.epochUsesTargetModel[i] = true;
                    matches++;
                }
            }
            if ( matches == 0 ) {
                throw new RuntimeException("Cannot find specified GLM substitution model (id: " + glmSubstitutionModel.getId() + ") in specified epoch model (id: " + epochBranchModel.getId() + ")");
            }
        } else {
            this.epochUsesTargetModel = new boolean[0];
        }

        this.nSites = nSites;
        this.nStates = glmSubstitutionModel.getFrequencyModel().getDataType().getStateCount();
        this.dim = nStates * (nStates - 1);

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

        this.proxyRates = new Parameter.Default(dim);
        this.proxy = new ComplexSubstitutionModel("internalGlmProxyForSubstitutionModelRandomEffectClassifier",
                glmSubstitutionModel.getDataType(),glmSubstitutionModel.getFrequencyModel(),proxyRates);
        this.markovJumps = new MarkovJumpsSubstitutionModel(proxy, MarkovJumpsType.COUNTS);
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

    private void makeProxyModel(int index, boolean includeRandomEffect) {
        double[] relativeRates = new double[dim];
        glmSubstitutionModel.setupRelativeRates(relativeRates);

        if (!includeRandomEffect) {
            double[] copiedParameterValues = glmSubstitutionModel.getGeneralizedLinearModel().getRandomEffect(0).getParameterValues();
//            double randomEffect = glmSubs.getGLM().getRandomEffect(0).getParameterValue(index);
            if (nullIsZero) {
                relativeRates[index] /= Math.exp(copiedParameterValues[index]);
            } else {
                relativeRates[index] = 0.0;
            }
        }

        for (int i = 0; i < dim; i++) {
            proxyRates.setParameterValue(i,relativeRates[i]);
        }

//        return new ComplexSubstitutionModel("internalGlmProxyForSubstitutionModelRandomEffectClassifier",
//                glmSubstitutionModel.getDataType(),glmSubstitutionModel.getFrequencyModel(),relativeRateDummyParameter);

    }

    private double getEpochContribution(NodeRef node) {
        BranchModel.Mapping map = epochBranchModel.getBranchModelMapping(node);
        int[] order = map.getOrder();
        double[] weights = map.getWeights();
        double duration = 0.0;
        for (int k = 0; k < order.length; k++) {
            if (epochUsesTargetModel[order[k]]) {
                duration += weights[k];
            }
        }
        return duration;
    }

    private double getTreeLengthInSubstitutions() {
        double length = 0.0;
        NodeRef root = tree.getRoot();
        for (int i = 0; i < tree.getNodeCount(); i++) {
            NodeRef node = tree.getNode(i);
            if ( node != root ) {
                if (usingEpochs) {
                    // This assumes branches do not have time-dependent rates that could interact with the epochs
                    length += getEpochContribution(node);
                } else {
                    double branchLength = tree.getNodeHeight(tree.getParent(node)) - tree.getNodeHeight(node);
                    length += branchLength * branchRates.getBranchRate(tree, node);
                }
            }
        }
        return length;
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
    
    private double getCountForRateCategory(int index, double time, boolean includeRandomEffect, boolean countAll) {
        makeProxyModel(index, includeRandomEffect);

        double[] register = new double[nStates * nStates];
        double[] jointCounts = new double[nStates * nStates];
        double[] conditionalCounts = new double[nStates * nStates];
        double[] transitionProbabilities = new double[nStates * nStates];

        if ( !countAll ) {
            int from = fromState[index];
            int to = toState[index];
            MarkovJumpsCore.fillRegistrationMatrix(register, from, to, nStates, 1.0);
        } else {
            MarkovJumpsCore.fillRegistrationMatrix(register, nStates);
        }
        markovJumps.setRegistration(register);

        markovJumps.computeJointStatMarkovJumps(time, jointCounts);
        markovJumps.computeCondStatMarkovJumps(time, conditionalCounts);
        proxy.getTransitionProbabilities(time, transitionProbabilities);

        double count = 0.0;
        double[] frequencies = glmSubstitutionModel.getFrequencyModel().getFrequencies();
        for (int i = 0; i < nStates; i++) {
            for (int j = 0; j < nStates; j++) {
//                count += jointCounts[i*nStates + j] * transitionProbabilities[i*nStates + j] * frequencies[i];
                count += jointCounts[i*nStates + j] * frequencies[i];
            }
        }

        return count;
//        return jointCounts[from * nStates + to];
    }
    
    private double getCount(int index, boolean includeRandomEffect, boolean countAll) {
        double duration = getTreeLengthInSubstitutions();
        double expectedCount = 0.0;
        if (usingRateVariation) {
            for (int i = 0; i < siteModel.getCategoryCount(); i++) {
                double thisRateCategory = siteModel.getRateForCategory(i);
                expectedCount += getCountForRateCategory(index, duration * thisRateCategory, includeRandomEffect, countAll)  * siteModel.getProportionForCategory(i);
            }
        } else {
            expectedCount += getCountForRateCategory(index, duration, includeRandomEffect, countAll);
        }

        return expectedCount * nSites;
    }

    private double getCountDifferences(int index) {
        return getCount(index,true, false) - getCount(index, false, false);
    }

    @Override
    public double getStatisticValue(int dim) {
//        double[] countAllSeparately = new double[this.dim];
//        double totCountAllSeparately = 0.0;
//        double totCountAll = getCount(-1, true, true);
//        for (int i = 0; i < this.dim; i++) {
//            countAllSeparately[i] = getCount(i,true, false);
//            totCountAllSeparately += countAllSeparately[i];
//        }
//
//        System.err.println("Total number of substitutions, counted separately: " + totCountAllSeparately + "\n");
//        System.err.println("Total number of substitutions: " + totCountAll + "\n");
//
//        double meanTreeLength = 0.0;
//        for (int i = 0; i < siteModel.getCategoryCount(); i++) {
//            meanTreeLength += siteModel.getProportionForCategory(i) * getTreeLengthInSubstitutions(siteModel.getRateForCategory(i));
//        }
//        System.err.println("meanTreeLength * nSites = " + meanTreeLength * nSites + "\n");

        return getDoubleResult(getCountDifferences(dim));
    }

}
