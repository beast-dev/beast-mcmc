/*
 * MarkovJumpsBeagleTreeLikelihood.java
 *
 * Copyright (C) 2002-2012 Alexei Drummond, Andrew Rambaut & Marc A. Suchard
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

package dr.app.beagle.evomodel.treelikelihood;

import dr.app.beagle.evomodel.branchmodel.BranchModel;
import dr.app.beagle.evomodel.sitemodel.BranchSubstitutionModel;
import dr.app.beagle.evomodel.sitemodel.SiteRateModel;
import dr.app.beagle.evomodel.substmodel.MarkovJumpsSubstitutionModel;
import dr.app.beagle.evomodel.substmodel.SubstitutionModel;
import dr.app.beagle.evomodel.substmodel.UniformizedSubstitutionModel;
import dr.evolution.alignment.PatternList;
import dr.evolution.datatype.DataType;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treelikelihood.TipStatesModel;
import dr.inference.loggers.LogColumn;
import dr.inference.loggers.NumberColumn;
import dr.inference.markovjumps.MarkovJumpsRegisterAcceptor;
import dr.inference.markovjumps.MarkovJumpsType;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;

import java.util.*;

/**
 * @author Marc Suchard
 * @author Vladimir Minin
 *         <p/>
 *         A base class for implementing Markov chain-induced counting processes (markovjumps) in BEAST using BEAGLE
 *         This work is supported by NSF grant 0856099
 *         <p/>
 *         Minin VN and Suchard MA (2008) Counting labeled transitions in continous-time Markov models of evolution.
 *         Journal of Mathematical Biology, 56, 391-412.
 */
public class MarkovJumpsBeagleTreeLikelihood extends AncestralStateBeagleTreeLikelihood
        implements MarkovJumpsRegisterAcceptor {

    public MarkovJumpsBeagleTreeLikelihood(PatternList patternList, TreeModel treeModel,
                                           BranchModel branchModel,
                                           SiteRateModel siteRateModel,
                                           BranchRateModel branchRateModel,
                                           TipStatesModel tipStatesModel,
                                           boolean useAmbiguities,
                                           PartialsRescalingScheme scalingScheme,
                                           Map<Set<String>, Parameter> partialsRestrictions,
                                           DataType dataType, String stateTag,
                                           boolean useMAP,
                                           boolean returnMarginalLikelihood,
                                           boolean useUniformization,
                                           boolean reportUnconditionedColumns,
                                           int nSimulants) {

        super(patternList, treeModel, branchModel, siteRateModel, branchRateModel, tipStatesModel, useAmbiguities,
                scalingScheme, partialsRestrictions, dataType, stateTag, useMAP, returnMarginalLikelihood);

        this.useUniformization = useUniformization;
        this.reportUnconditionedColumns = reportUnconditionedColumns;
        this.nSimulants = nSimulants;

        markovjumps = new ArrayList<MarkovJumpsSubstitutionModel>();
        branchModelNumber = new ArrayList<Integer>();
        registerParameter = new ArrayList<Parameter>();
        jumpTag = new ArrayList<String>();
        expectedJumps = new ArrayList<double[][]>();
//        storedExpectedJumps = new ArrayList<double[][]>();

        tmpProbabilities = new double[stateCount * stateCount * categoryCount];
        condJumps = new double[categoryCount][stateCount * stateCount];
    }

    public void addRegister(Parameter addRegisterParameter,
                            MarkovJumpsType type,
                            boolean scaleByTime) {

        if ((type == MarkovJumpsType.COUNTS &&
                addRegisterParameter.getDimension() != stateCount * stateCount) ||
                (type == MarkovJumpsType.REWARDS &&
                        addRegisterParameter.getDimension() != stateCount)
                ) {
            throw new RuntimeException("Register parameter of wrong dimension");
        }
        addVariable(addRegisterParameter);
        final String tag = addRegisterParameter.getId();

        for (int i = 0; i < substitutionModelDelegate.getSubstitutionModelCount(); ++i) {

            registerParameter.add(addRegisterParameter);
            MarkovJumpsSubstitutionModel mjModel;
            SubstitutionModel substitutionModel = substitutionModelDelegate.getSubstitutionModel(i);

            if (useUniformization) {
                mjModel = new UniformizedSubstitutionModel(substitutionModel, type, nSimulants);
            } else {
                if (type == MarkovJumpsType.HISTORY) {
                    throw new RuntimeException("Can only report complete history using uniformization");
                }
                mjModel = new MarkovJumpsSubstitutionModel(substitutionModel, type);
            }
            markovjumps.add(mjModel);
            branchModelNumber.add(i);
            addModel(mjModel);
            setupRegistration(numRegisters);

            String traitName;

            if (substitutionModelDelegate.getSubstitutionModelCount() == 1) {
                traitName = tag;
            } else {
                traitName = tag + i;
            }

            jumpTag.add(traitName);
            
            expectedJumps.add(new double[treeModel.getNodeCount()][patternCount]);
//        storedExpectedJumps.add(new double[treeModel.getNodeCount()][patternCount]);

            boolean[] oldScaleByTime = this.scaleByTime;
            int oldScaleByTimeLength = (oldScaleByTime == null ? 0 : oldScaleByTime.length);
            this.scaleByTime = new boolean[oldScaleByTimeLength + 1];
            if (oldScaleByTimeLength > 0) {
                System.arraycopy(oldScaleByTime, 0, this.scaleByTime, 0, oldScaleByTimeLength);
            }
            this.scaleByTime[oldScaleByTimeLength] = scaleByTime;

            if (type != MarkovJumpsType.HISTORY) {
                treeTraits.addTrait(traitName, new TreeTrait.DA() {

                    final int registerNumber = numRegisters;

                    public String getTraitName() {
                        return tag;
                    }

                    public Intent getIntent() {
                        return Intent.BRANCH;
                    }

                    public double[] getTrait(Tree tree, NodeRef node) {
                        return getMarkovJumpsForNodeAndRegister(tree, node, registerNumber);
                    }
                });
            } else {
                if (histories == null) {
                    histories = new String[treeModel.getNodeCount()][patternCount];
                } else {
                    throw new RuntimeException("Only one complete history per markovJumpTreeLikelihood is allowed");
                }
                if (nSimulants > 1) {
                    throw new RuntimeException("Only one simulant allowed when saving complete history");
                }

                // Add total number of changes over tree trait
                TreeTrait da = new TreeTrait.DA() {

                    final int registerNumber = numRegisters;

                    public String getTraitName() {
                        return tag;
                    }

                    public Intent getIntent() {
                        return Intent.BRANCH;
                    }

                    public double[] getTrait(Tree tree, NodeRef node) {
                        return getMarkovJumpsForNodeAndRegister(tree, node, registerNumber);
                    }
                };

                treeTraits.addTrait(addRegisterParameter.getId(), new TreeTrait.SumOverTreeDA(da));

                historyRegisterNumber = numRegisters; // Record the complete history for this register
                ((UniformizedSubstitutionModel) mjModel).setSaveCompleteHistory(true);

                for (int site = 0; site < patternCount; ++site) {

                    final String anonName = (patternCount == 1) ? HISTORY : HISTORY + "_" + (site + 1);
                    final int anonSite = site;

                    treeTraits.addTrait(anonName, new TreeTrait.S() {

                        public String getTraitName() {
                            return anonName;
                        }

                        public Intent getIntent() {
                            return Intent.BRANCH;
                        }

                        public String getTrait(Tree tree, NodeRef node) {
                            return getHistoryForNode(tree, node, anonSite);
                        }

                        public boolean getLoggable() {
                            return logHistory;
                        }
                    });
                }
            }

            numRegisters++;

        } // End of loop over branch models
    }

    public void setLogHistories(boolean in) {
        logHistory = in;
    }

//    public double[] getRewardsForNodeAndPattern(Tree tree, NodeRef node, int pattern) {
//        double[] rtn = new double[numRegisters];
//        for (int r = 0; r < numRegisters; r++) {
//            double[] mjs = getMarkovJumpsForNodeAndRegister(tree, node, r);
//            rtn[r] = mjs[pattern];
//        }
//        return rtn;
//    }


    public double[] getMarkovJumpsForNodeAndRegister(Tree tree, NodeRef node, int whichRegister) {
        return getMarkovJumpsForRegister(tree, whichRegister)[node.getNumber()];
    }

//    public double[][] getMarkovJumpsForNode(Tree tree, NodeRef node) {
//        double[][] rtn = new double[numRegisters][];
//        for (int r = 0; r < numRegisters; r++) {
//            rtn[r] = getMarkovJumpsForNodeAndRegister(tree, node, r);
//        }
//        return rtn;
//    }

    private void refresh(Tree tree) {
        if (tree != treeModel) {
            throw new RuntimeException("Must call with internal tree");
        }

        if (!likelihoodKnown) {
            calculateLogLikelihood();
            likelihoodKnown = true;
        }
        
        if (!areStatesRedrawn) {
            redrawAncestralStates();
        }
    }

    public double[][] getMarkovJumpsForRegister(Tree tree, int whichRegister) {
        refresh(tree);
        return expectedJumps.get(whichRegister);
    }

    public String getHistoryForNode(Tree tree, NodeRef node, int site) {
        return getHistory(tree)[node.getNumber()][site];
    }

    public String[][] getHistory(Tree tree) {
        refresh(tree);
        return histories;
    }

//    private static String formattedValue(double[] values) {
//        double total = 0;
//        for (double summant : values) {
//            total += summant;
//        }
//        return Double.toString(total); // Currently return the sum across sites
//    }

    private void setupRegistration(int whichRegistration) {

        double[] registration = registerParameter.get(whichRegistration).getParameterValues();
        markovjumps.get(whichRegistration).setRegistration(registration);
        areStatesRedrawn = false;
    }

    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        for (int r = 0; r < numRegisters; r++) {
            if (variable == registerParameter.get(r)) {
                setupRegistration(r);
                return;
            }
        }
        super.handleVariableChangedEvent(variable, index, type);
    }

    protected void hookCalculation(Tree tree, NodeRef parentNode, NodeRef childNode,
                                   int[] parentStates, int[] childStates,
                                   double[] inProbabilities, int[] rateCategory) {

        final int childNum = childNode.getNumber();

        double[] probabilities = inProbabilities;
        if (probabilities == null) { // Leaf will call this hook with a null
            getMatrix(childNum, tmpProbabilities);
            probabilities = tmpProbabilities;
        }

        final double branchRate = branchRateModel.getBranchRate(tree, childNode);
        final double parentTime = tree.getNodeHeight(parentNode);
        final double childTime = tree.getNodeHeight(childNode);
        final double substTime = parentTime - childTime;

        for (int r = 0; r < markovjumps.size(); r++) {
            MarkovJumpsSubstitutionModel thisMarkovJumps = markovjumps.get(r);

            final int modelNumberFromrRegistry = branchModelNumber.get(r);
//            int dummy = 0;
//            final int modelNumberFromTree = branchSubstitutionModel.getBranchIndex(tree, childNode, dummy);
            // @todo AR - not sure about this - if this is an epoch this is just going to get the most
            // @todo tipward model for the branch. I think this was what was happening before (in comment,
            // @todo above).
            BranchModel.Mapping mapping = branchModel.getBranchModelMapping(childNode);

            if (modelNumberFromrRegistry == mapping.getOrder()[0]) {
                if (useUniformization) {
                    computeSampledMarkovJumpsForBranch(((UniformizedSubstitutionModel) thisMarkovJumps), substTime,
                            branchRate, childNum, parentStates, childStates, parentTime, childTime,probabilities, scaleByTime[r],
                            expectedJumps.get(r), rateCategory, r == historyRegisterNumber);
                } else {
                    computeIntegratedMarkovJumpsForBranch(thisMarkovJumps, substTime, branchRate, childNum, parentStates,
                            childStates, probabilities, condJumps, scaleByTime[r], expectedJumps.get(r), rateCategory);
                }
            } else {
                // Fill with zeros
                double[] result = expectedJumps.get(r)[childNum];
                Arrays.fill(result, 0.0);
            }
        }
    }

    private void computeSampledMarkovJumpsForBranch(UniformizedSubstitutionModel thisMarkovJumps,
                                                    double substTime,
                                                    double branchRate,
                                                    int childNum,
                                                    int[] parentStates,
                                                    int[] childStates,
                                                    double parentTime,
                                                    double childTime,
                                                    double[] probabilities,
                                                    boolean scaleByTime,
                                                    double[][] thisExpectedJumps,
                                                    int[] rateCategory,
                                                    boolean saveHistory) {


        // Fill condJumps with sampled values for this branch for each site
        for (int j = 0; j < patternCount; j++) {
            final int category = rateCategory == null ? 0 : rateCategory[j];
            final double categoryRate = siteRateModel.getRateForCategory(category);
            final int matrixIndex = category * stateCount * stateCount;
            double value = thisMarkovJumps.computeCondStatMarkovJumps(
                    parentStates[j],
                    childStates[j],
                    substTime * branchRate * categoryRate,
                    probabilities[matrixIndex + parentStates[j] * stateCount + childStates[j]]
            );
            if (scaleByTime) {
                value /= branchRate * categoryRate;
            }
            thisExpectedJumps[childNum][j] = value;
            if (saveHistory) {
                histories[childNum][j] = thisMarkovJumps.getCompleteHistory(parentTime, childTime);
            }
        }
    }

    private void computeIntegratedMarkovJumpsForBranch(MarkovJumpsSubstitutionModel thisMarkovJumps,
                                                       double substTime,
                                                       double branchRate,
                                                       int childNum,
                                                       int[] parentStates,
                                                       int[] childStates,
                                                       double[] probabilities,
                                                       double[][] condJumps,
                                                       boolean scaleByTime,
                                                       double[][] thisExpectedJumps,
                                                       int[] rateCategory) {

        // Fill condJumps with conditional mean values for this branch
        for (int i = 0; i < categoryCount; i++) {

            double rate = siteRateModel.getRateForCategory(i);
            if (rate > 0) {

                if (categoryCount == 1) {
                    thisMarkovJumps.computeCondStatMarkovJumps(
                            substTime * branchRate * rate,
                            probabilities,
                            condJumps[i]);
                } else {

                    System.arraycopy(probabilities, i * stateCount * stateCount, tmpProbabilities, 0,
                            stateCount * stateCount);

                    thisMarkovJumps.computeCondStatMarkovJumps(
                            substTime * branchRate * rate,
                            tmpProbabilities,
                            condJumps[i]);
                }

                if (scaleByTime) {
                    double scalar = branchRate * rate;
                    for (int j = 0; j < condJumps[i].length; j++) {
                        condJumps[i][j] /= scalar;
                    }
                }
            } else {
                Arrays.fill(condJumps[i], 0.0);
                if (thisMarkovJumps.getType() == MarkovJumpsType.REWARDS && scaleByTime) {
                    for (int j = 0; j < stateCount; j++) {
                        condJumps[i][j * stateCount + j] = substTime;
                    }
                }
            }
        }

        for (int j = 0; j < patternCount; j++) { // Pick out values given parent and child states
            int category = rateCategory == null ? 0 : rateCategory[j];
            thisExpectedJumps[childNum][j] = condJumps[category][parentStates[j] * stateCount + childStates[j]];
        }
    }

//    public void storeState() {
//
//        super.storeState();
//
////        if (areStatesRedrawn) {
//            for (int i = 0; i < expectedJumps.size(); i++) {
//                double[][] thisExpectedJumps = expectedJumps.get(i);
//                double[][] storedThisExpectedJumps = storedExpectedJumps.get(i);
//                for (int j = 0; j < thisExpectedJumps.length; j++) {
//                    System.arraycopy(thisExpectedJumps[j], 0, storedThisExpectedJumps[j], 0,
//                            thisExpectedJumps[j].length);
//                }
//            }
////        }
//    }
//
//    public void restoreState() {
//
//        super.restoreState();
//
//        List<double[][]> tmp = expectedJumps;
//        expectedJumps = storedExpectedJumps;
//        storedExpectedJumps = tmp;
//        areStatesRedrawn = false;
//
//    }

    public LogColumn[] getColumns() {

        int nColumns = patternCount * numRegisters;
        if (reportUnconditionedColumns) {
            if (categoryCount == 1) {
                nColumns += numRegisters;
            } else {
                nColumns *= 2;
            }
        }

        int index = 0;
        LogColumn[] allColumns = new LogColumn[nColumns];
        for (int r = 0; r < numRegisters; r++) {
            for (int j = 0; j < patternCount; j++) {
                allColumns[index++] = new ConditionedCountColumn(jumpTag.get(r), r, j);
                if (reportUnconditionedColumns) {
                    if (categoryCount > 1) {
                        allColumns[index++] = new UnconditionedCountColumn(jumpTag.get(r), r, j, rateCategory);
                    }
                }
            }
            if (reportUnconditionedColumns) {
                if (categoryCount == 1) {
                    allColumns[index++] = new UnconditionedCountColumn(jumpTag.get(r), r);
                }
            }
        }
        return allColumns;
    }

    protected abstract class CountColumn extends NumberColumn {
        protected int indexRegistration;
        protected int indexSite;

        public CountColumn(String label, int r, int j) {
            super(label + (j >= 0 ? "[" + (j + 1) + "]" : ""));
            indexRegistration = r;
            indexSite = j;
        }

        public abstract double getDoubleValue();
    }

    protected class ConditionedCountColumn extends CountColumn {

        public ConditionedCountColumn(String label, int r, int j) {
            super("c_" + label, r, j);
        }

        public double getDoubleValue() {
            double total = 0;
            double[][] values = getMarkovJumpsForRegister(treeModel, indexRegistration);
            for (int i = 0; i < treeModel.getNodeCount(); i++) {
                total += values[i][indexSite];
            }
            return total;
        }
    }

    protected class UnconditionedCountColumn extends CountColumn {
        int[] rateCategory;

        public UnconditionedCountColumn(String label, int r, int j, int[] rateCategory) {
            super("u_" + label, r, j);
            this.rateCategory = rateCategory;
        }

        public UnconditionedCountColumn(String label, int r) {
            this(label, r, -1, null);
        }

        public double getDoubleValue() {
            double value = markovjumps.get(indexRegistration).getMarginalRate() * getExpectedTreeLength();
            if (rateCategory != null) {
                value *= siteRateModel.getRateForCategory(rateCategory[indexSite]);
            }
            return value;
        }

        private double getExpectedTreeLength() {
            double expectedTreeLength = 0;
            for (int i = 0; i < treeModel.getNodeCount(); i++) {
                NodeRef node = treeModel.getNode(i);
                if (!treeModel.isRoot(node)) {
                    expectedTreeLength += branchRateModel.getBranchRate(treeModel, node)
                            * treeModel.getBranchLength(node);
                }
            }
            return expectedTreeLength;
        }
    }

    public static final String HISTORY = "history";
    public static final String TOTAL_COUNTS = "allTransitions";

    private List<MarkovJumpsSubstitutionModel> markovjumps;
    private List<Integer> branchModelNumber;
    private List<Parameter> registerParameter;
    private List<String> jumpTag;
    private List<double[][]> expectedJumps;
//    private List<double[][]> storedExpectedJumps;
    private boolean logHistory = false;
    private String[][] histories = null;
    private boolean[] scaleByTime;
    private double[] tmpProbabilities;
    private double[][] condJumps;
    private int numRegisters;
    private int historyRegisterNumber = -1;
    private final boolean useUniformization;
    private final int nSimulants;
    private final boolean reportUnconditionedColumns;
}
