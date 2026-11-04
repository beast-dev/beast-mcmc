/*
 * EpochBranchModel.java
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

package dr.evomodel.branchmodel;

import dr.evolution.datatype.DataType;
import dr.evolution.tree.NodeRef;
import dr.evomodel.branchratemodel.ArbitraryBranchRates;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.substmodel.*;
import dr.evomodel.tree.TreeModel;
import dr.inference.markovjumps.SericolaSeriesMarkovReward;
import dr.inference.model.*;
import dr.util.Author;
import dr.util.Citable;
import dr.util.Citation;
import dr.xml.Reportable;

import java.util.*;

/**
 * @author Filippo Monti
 * @author Marc A. Suchard
 */

public class RewardsAwareBranchModel extends AbstractModel implements TransitionMatrixProviderBranchModel, Citable, Reportable {

    public static final String REWARDS_AWARE_BRANCH_MODEL = "RewardsAwareBranchModel";
    private final Parameter rewardRates;
    private final SubstitutionModel underlyingSubstitutionModel;
    private final TreeModel tree;
    private final ArbitraryBranchRates branchRateModel;

    private final int nstates;
    private double[] unsortedQ;
    private double[][] unsortedW;

    private double[] sortedRewardRates;

    private double[] Q;
    private double[][] W;

    private int[] perm;
    private List<SubstitutionModel> substitutionModels;
    private boolean knownTransitionMatrices;
    private boolean knownSortedQ;
    private boolean knownSortedRewardRates;

    public RewardsAwareBranchModel(TreeModel tree,
                                   SubstitutionModel underlyingSubstitutionModel, //TODO askMarc why LogRateSubModel does not work
                                   Parameter rewardRates,
                                   ArbitraryBranchRates branchRateModel) {

        super(REWARDS_AWARE_BRANCH_MODEL);

        this.underlyingSubstitutionModel = underlyingSubstitutionModel;
        this.rewardRates = rewardRates;
        this.branchRateModel = branchRateModel;

        this.tree = tree;
        if (underlyingSubstitutionModel == null) {
            throw new IllegalArgumentException("RewardsAwareBranchModel must be provided " +
                    "with an underlying substitution model");
        }

        this.nstates = underlyingSubstitutionModel.getDataType().getStateCount();
        this.sortedRewardRates = new double[nstates];

        this.unsortedQ = new double[nstates * nstates];
        this.Q = new double[nstates * nstates];

        this.W = new double[tree.getNodeCount()][nstates * nstates];
        this.unsortedW = new double[tree.getNodeCount()][nstates * nstates];

        this.perm = new int[nstates];

        this.knownSortedRewardRates = false;
        this.knownSortedQ = false;
        this.knownTransitionMatrices = false;

        addModel(underlyingSubstitutionModel);
        addModel(branchRateModel);
        addVariable(rewardRates);
    }

    public FrequencyModel getRootFrequencyModel() {return underlyingSubstitutionModel.getFrequencyModel();}

    @Override
    public SubstitutionModel getRootSubstitutionModel() {return underlyingSubstitutionModel;}

    boolean DUMMYTESTING = false;
    boolean DEBUG = false;

    public double[] getTransitionMatrix(int i) {
        if (DUMMYTESTING) { //TODO this is only for testing
            NodeRef node = tree.getNode(i);
            double branchLength = tree.getBranchLength(node);
            getRootSubstitutionModel().getTransitionProbabilities(branchLength, W[i]);
        } else {
            computeTransitionMatrices();
        }
        return W[i];
    }

    @Override
    public double[] getTransitionMatrix(NodeRef branch) {
        return getTransitionMatrix(branch.getNumber()); // TODO this should be direct
    }

    private void computeTransitionMatrices() {
        if (DEBUG) {System.out.println("computeTransitionMatrices");}
        if (!knownTransitionMatrices) {
            sortRewardRates();
            sortQ();
            SericolaSeriesMarkovReward sericolaSeriesMarkovReward =
                    new SericolaSeriesMarkovReward(Q, sortedRewardRates, nstates);

//            final int branchCount = tree.getNodeCount() - 1;
//            double[] branchRates = new double[branchCount];
//            double[] branchLengths = new double[branchCount];

//            int k = 0;
            for (int i = 0; i < tree.getNodeCount(); i++) {
                NodeRef node = tree.getNode(i);
                if (tree.isRoot(node)) continue;
                final double branchRate = branchRateModel.getBranchRate(tree, node);
                final double branchLength = tree.getBranchLength(node);
//                branchRates[k] = branchRateModel.getBranchRate(tree, node);
//                branchLengths[k] = tree.getBranchLength(node);
//                nodeIndicesW[node.getNumber()] = k;
                unsortedW[i] = sericolaSeriesMarkovReward.computePdf(branchRate, branchLength);
//                k++;
            }
//            unsortedW = sericolaSeriesMarkovReward.computePdf(branchRates, branchLengths, true); // TODO update using this call directly
            sortW();
            knownTransitionMatrices = true;
        }
    }

    private void sortRewardRates() { // compute permutation and sort reward rates
        if (!knownSortedRewardRates) {
            System.out.println("Sorting reward rates");
            final double[] rewardVals = rewardRates.getParameterValues();

            Integer[] permObj = new Integer[nstates];
            for (int i = 0; i < nstates; i++) permObj[i] = i;

            Arrays.sort(permObj, Comparator.comparingDouble(i -> rewardVals[i]));

            for (int i = 0; i < nstates; i++) {
                perm[i] = permObj[i];
                sortedRewardRates[i] = rewardVals[perm[i]];
            }
            knownSortedQ = false;
            knownSortedRewardRates = true;
        }
    }

    private void sortQ() {
        if (!knownSortedQ) {
            underlyingSubstitutionModel.getInfinitesimalMatrix(unsortedQ);
            //  reorder Q using perm (column-major)
            for (int j = 0; j < nstates; j++) {
                int pj = perm[j];
                for (int i = 0; i < nstates; i++) {
                    Q[j * nstates + i] = unsortedQ[pj * nstates + perm[i]];
                }
            }
            knownSortedQ = true;
        }
    }

    private void sortW() {
        // Build inverse permutation
        int[] invPerm = new int[nstates];
        for (int i = 0; i < nstates; i++) {
            invPerm[perm[i]] = i;
        }

        // Apply inverse permutation to recover original order
        for (int k = 0; k < tree.getNodeCount(); k++) {
            if (tree.isRoot(tree.getNode(k))) continue;
//            int nodeIndex = nodeIndicesW[k];
            double[] wb = unsortedW[k];
            double[] swb = W[k];
            for (int j = 0; j < nstates; j++) {
                final int ij = invPerm[j] * nstates;
                final int jj = j * nstates;
                for (int i = 0; i < nstates; i++) {
                    swb[jj + i] = wb[ij + invPerm[i]];
                }
            }
        }
    }

    @Override
    public Mapping getBranchModelMapping(NodeRef node) {
        // TODO this is redundant: it should map directly to the finite time transition matrix
        final double[] weights = new double[1];
        weights[0] = 1.0;
        final int[] order = new int[1];
        order[0] = node.getNumber();
//        order[0] = nodeIndicesW[]

        return new Mapping() {
            @Override
            public int[] getOrder() {
                return order;
            }

            @Override
            public double[] getWeights() {
                return weights;
            }
        };
    }

    @Override
    public boolean requiresMatrixConvolution() {return false;}

    public TreeModel getTree() { return tree; }
    public Parameter getRewardRates() { return rewardRates; }
    public BranchRateModel getRateBranchModel() { return branchRateModel; }


    boolean ignoreModelChangedEvent = false;

    protected void handleModelChangedEvent(Model model, Object object, int index) {
        if (ignoreModelChangedEvent) {return;}
        System.out.println("Model changed event");
        if (model == underlyingSubstitutionModel) {
            knownSortedQ = false;
            knownTransitionMatrices = false;
            fireModelChanged();
        } else if (model == branchRateModel) { // the branch rate model provides the tree and total rewards
            knownSortedQ = false;
            knownTransitionMatrices = false;
            fireModelChanged();
        } else {
            throw new IllegalArgumentException("Unknown model: " + model);
        }
    }

    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        if (variable == rewardRates) {
            System.out.println("Reward rates changed");
            knownSortedRewardRates = false;
            knownTransitionMatrices = false;
            //TODO this is not firing a change in the treeDataLikelihood that build again the transition matrices
            fireModelChanged();
        } else {
            throw new IllegalArgumentException("Unknown variable: " + variable);
        }
    }

    double[][] storedW;
    protected void storeState() {
        storedW = Arrays.copyOf(W, W.length);
    }

    protected void restoreState() {
        // TODO save Q and sortedRewardRates?
        double[][] tempW = storedW;
        storedW = W;
        W = tempW;
    }

    protected void acceptState() {
    }

    @Override
    public Citation.Category getCategory() {
        return Citation.Category.SUBSTITUTION_MODELS;
    }

    @Override
    public String getDescription() {
        return "Rewards Aware Branch model";
    }

    @Override
    public List<Citation> getCitations() {
        return Collections.singletonList(
                new Citation(new Author[]{new Author("F", "Monti"),
                        new Author("MA", "Suchard")},
                        "Dependencies between CTMCs",
                        2025,
                        "TOBE",
                        1,
                        1,
                        1,
                        Citation.Status.IN_PRESS));
    }

    @Override
    public String getReport() {
        if (!DUMMYTESTING) {
            computeTransitionMatrices();
        }
        StringBuilder sb = new StringBuilder();
        sb.append("W matrix: ");
        for (int i=0; i < tree.getNodeCount() - 1; i++) {
//            sb.append("Branch ").append(i).append(": ");
            for (double val : W[i]) {
                sb.append(val).append(" ");
            }
//            sb.append(Arrays.toString(W[i]));
//                    .append("\n");
        }
        return sb.toString();
    }

    @Override
    public List<SubstitutionModel> getSubstitutionModels() {
        System.out.println("getSubstitutionModels");
        buildSubstitutionModels(); // TODO this is only for compatibility with SubstitutionModelDelegate; It is called only once
        return substitutionModels;
    }

    protected void buildSubstitutionModels() {

//        computeTransitionMatrices();

        substitutionModels = new ArrayList<>();

        for (int i = 0; i < tree.getNodeCount(); i++) {
            NodeRef node = tree.getNode(i);
            if (tree.isRoot(node)) continue;
            ignoreModelChangedEvent = true; // todo this can be deleted since this method is called only once
            SubstitutionModel substitutionModel = new TransitionMatrixProvider(
                    "RewardsAwareSubstitutionModel",
                    underlyingSubstitutionModel.getDataType(),
                    underlyingSubstitutionModel.getFrequencyModel(),
//                    W[nodeIndicesW[node.getNumber()]]
                    W[node.getNumber()]
            );
            ignoreModelChangedEvent = false;
            substitutionModels.add(substitutionModel);
        }
    }

    class TransitionMatrixProvider extends ComplexSubstitutionModel { // TODO this should be deprecated
        private double[] transitionMatrix;

        public TransitionMatrixProvider(String name, DataType dataType,
                                        FrequencyModel freqModel, double[] transitionMatrix) {
            super(name, dataType, freqModel, null);
            this.transitionMatrix = transitionMatrix;
        }

        @Override
        public void getTransitionProbabilities(double distance, double[] matrix) {
            System.arraycopy(transitionMatrix, 0, matrix, 0, transitionMatrix.length);
        }

        protected void handleModelChangedEvent(Model model, Object object, int index) {}
        @Override
        protected void frequenciesChanged() {}
        @Override
        protected void ratesChanged() {}
        @Override
        protected void setupRelativeRates(double[] rates) {}

    }
}

