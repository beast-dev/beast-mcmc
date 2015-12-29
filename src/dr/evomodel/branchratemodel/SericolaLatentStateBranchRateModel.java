/*
 * SericolaLatentStateBranchRateModel.java
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

package dr.evomodel.branchratemodel;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.tree.TreeParameterModel;
import dr.inference.markovjumps.MarkovReward;
import dr.inference.markovjumps.TwoStateOccupancyMarkovReward;
import dr.inference.model.AbstractModelLikelihood;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;

import java.util.ArrayList;
import java.util.List;

/**
 * SericolaLatentStateBranchRateModel
 *
 * @author Andrew Rambaut
 * @author Marc Suchard
 * @version $Id$
 *          <p/>
 *          $HeadURL$
 *          <p/>
 *          $LastChangedBy$
 *          $LastChangedDate$
 *          $LastChangedRevision$
 */
public class SericolaLatentStateBranchRateModel extends AbstractModelLikelihood implements BranchRateModel {

    public static final String LATENT_STATE_BRANCH_RATE_MODEL = "latentStateBranchRateModel";

    public static final boolean USE_CACHING = true;
    // seed 666, caching off: 204.69 seconds for 20000 states
    // state 20000	-5510.2520
    // 85.7%  5202  +     6    dr.inference.markovjumps.SericolaSeriesMarkovReward.accumulatePdf

    // seed 666, caching on: 119.43 seconds for 20000 states
    // state 20000	-5510.2520
    // 83.4%  3156  +     4    dr.inference.markovjumps.SericolaSeriesMarkovReward.accumulatePdf

    private final TreeModel tree;
    private final BranchRateModel nonLatentRateModel;
    private final Parameter latentTransitionRateParameter;
    private final Parameter latentTransitionFrequencyParameter;
    private final TreeParameterModel latentStateProportions;
    private final Parameter latentStateProportionParameter;
    private final CountableBranchCategoryProvider branchCategoryProvider;

    private MarkovReward series;
    private MarkovReward storedSeries;
    private boolean likelihoodKnown = false;
    private boolean storedLikelihoodKnown;
    private double logLikelihood;
    private double storedLogLikelihood;

    private double[] branchLikelihoods;
    private double[] storedbranchLikelihoods;

    private boolean[] updateBranch;
    private boolean[] storedUpdateBranch;

    private boolean[] updateCategory;
    private boolean[] storedUpdateCategory;

    public SericolaLatentStateBranchRateModel(String name,
                                              TreeModel treeModel,
                                              BranchRateModel nonLatentRateModel,
                                              Parameter latentTransitionRateParameter,
                                              Parameter latentTransitionFrequencyParameter,
                                              Parameter latentStateProportionParameter,
                                              CountableBranchCategoryProvider branchCategoryProvider) {
        super(name);

        this.tree = treeModel;
        addModel(tree);

        this.nonLatentRateModel = nonLatentRateModel;
        addModel(nonLatentRateModel);

        this.latentTransitionRateParameter = latentTransitionRateParameter;
        addVariable(latentTransitionRateParameter);

        this.latentTransitionFrequencyParameter = latentTransitionFrequencyParameter;
        addVariable(latentTransitionFrequencyParameter);

        if (branchCategoryProvider ==  null) {
            this.latentStateProportions = new TreeParameterModel(tree, latentStateProportionParameter, false, Intent.BRANCH);
            addModel(latentStateProportions);

            this.latentStateProportionParameter = null;
            this.branchCategoryProvider = null;
        } else {
            this.latentStateProportions = null;
            this.branchCategoryProvider = branchCategoryProvider;
            this.latentStateProportionParameter = latentStateProportionParameter;
            this.latentStateProportionParameter.setDimension(branchCategoryProvider.getCategoryCount());

            if (USE_CACHING) {
                updateCategory = new boolean[branchCategoryProvider.getCategoryCount()];
                storedUpdateCategory = new boolean[branchCategoryProvider.getCategoryCount()];
                setUpdateAllCategories();
            }

            addVariable(latentStateProportionParameter);
        }

        branchLikelihoods = new double[tree.getNodeCount()];
        if (USE_CACHING) {
            updateBranch = new boolean[tree.getNodeCount()];
            storedUpdateBranch = new boolean[tree.getNodeCount()];
            storedbranchLikelihoods = new double[tree.getNodeCount()];

            setUpdateAllBranches();
        }
    }

    public SericolaLatentStateBranchRateModel(Parameter rate, Parameter prop) {
        super(LATENT_STATE_BRANCH_RATE_MODEL);
        tree = null;
        nonLatentRateModel = null;
        latentTransitionRateParameter = rate;
        latentTransitionFrequencyParameter = prop;
        latentStateProportions = null;
        this.latentStateProportionParameter = null;
        this.branchCategoryProvider = null;
    }

    private double[] createLatentInfinitesimalMatrix() {
        final double rate = latentTransitionRateParameter.getParameterValue(0);
        final double prop = latentTransitionFrequencyParameter.getParameterValue(0);

        double[] mat = new double[]{
                -rate * prop, rate * prop,
                rate * (1.0 - prop), -rate * (1.0 - prop)
        };
        return mat;
    }

    private static double[] createReward() {
        return new double[]{0.0, 1.0};
    }

    private MarkovReward createSeries() {
//        MarkovReward series = new SericolaSeriesMarkovReward(createLatentInfinitesimalMatrix(),
//                createReward(), 2);
        MarkovReward series = new TwoStateOccupancyMarkovReward(createLatentInfinitesimalMatrix());
        return series;
    }

    @Override
    public double getBranchRate(Tree tree, NodeRef node) {
        double nonLatentRate = nonLatentRateModel.getBranchRate(tree, node);

        double latentProportion = getLatentProportion(tree, node);

        return calculateBranchRate(nonLatentRate, latentProportion);
    }

    public double getLatentProportion(Tree tree, NodeRef node) {

        if (latentStateProportions != null) {
            return latentStateProportions.getNodeValue(tree, node);
        } else {
            return latentStateProportionParameter.getParameterValue(branchCategoryProvider.getBranchCategory(tree, node));
        }
    }

    private double calculateBranchRate(double nonLatentRate, double latentProportion) {
        return nonLatentRate * (1.0 - latentProportion);
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        if (model == tree) {
            likelihoodKnown = false; // node heights change elapsed times on branches, TODO could cache

            if (index == -1) {
                setUpdateAllBranches();
            } else {
                setUpdateBranch(index);
            }

        } else if (model == nonLatentRateModel) {
            // rates will change but the latent proportions haven't so the density is unchanged
        } else if (model == latentStateProportions) {
            likelihoodKnown = false; // argument of density has changed

            if (index == -1) {
                setUpdateAllBranches();
            } else {
                setUpdateBranch(index);
            }
        }
        fireModelChanged();
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        if (variable == latentTransitionFrequencyParameter || variable == latentTransitionRateParameter) {
            // series computations have changed
            series = null;
            setUpdateAllBranches();
            likelihoodKnown = false;
        } else if (variable == latentStateProportionParameter) {
            if (index == -1) {
                setUpdateAllBranches();
            } else {
                setUpdateBranchCategory(index);
            }
            likelihoodKnown = false;
            fireModelChanged();
        }
    }

    private void setUpdateBranch(int nodeNumber) {
        if (USE_CACHING) {
            updateBranch[nodeNumber] = true;
        }
    }

    private void setUpdateAllBranches() {
        if (USE_CACHING) {
            for (int i = 0; i < updateBranch.length; i++) {
                updateBranch[i] = true;
            }
        }
    }

    private void clearUpdateAllBranches() {
        if (USE_CACHING) {
            for (int i = 0; i < updateBranch.length; i++) {
                updateBranch[i] = false;
            }
        }
    }

    private void setUpdateBranchCategory(int category) {
        if (USE_CACHING) {
            updateCategory[category] = true;

        }
    }

    private void setUpdateAllCategories() {
        if (USE_CACHING) {
            for (int i = 0; i < updateCategory.length; i++) {
                updateCategory[i] = true;
            }

        }
    }

    private void clearAllCategories() {
        if (USE_CACHING && updateCategory != null) {
            for (int i = 0; i < updateCategory.length; i++) {
                updateCategory[i] = false;
            }

        }
    }

    @Override
    protected void storeState() {
        storedSeries = series;
        storedLogLikelihood = logLikelihood;
        storedLikelihoodKnown = likelihoodKnown;
        if (USE_CACHING) {
            System.arraycopy(branchLikelihoods, 0, storedbranchLikelihoods, 0, branchLikelihoods.length);
            System.arraycopy(updateBranch, 0, storedUpdateBranch, 0, updateBranch.length);

            if (updateCategory != null) {
                System.arraycopy(updateCategory, 0, storedUpdateCategory, 0, updateCategory.length);
            }
        }
    }

    @Override
    protected void restoreState() {
        series = storedSeries;
        logLikelihood = storedLogLikelihood;
        likelihoodKnown = storedLikelihoodKnown;

        if (USE_CACHING) {
            double[] tmp = branchLikelihoods;
            branchLikelihoods = storedbranchLikelihoods;
            storedbranchLikelihoods = tmp;

            boolean[] tmp2 = updateBranch;
            updateBranch = storedUpdateBranch;
            storedUpdateBranch = tmp2;

            boolean[] tmp3 = updateCategory;
            updateCategory = storedUpdateCategory;
            storedUpdateCategory = tmp3;
        }
    }

    @Override
    protected void acceptState() {

    }

    @Override
    public Model getModel() {
        return this;
    }

    @Override
    public double getLogLikelihood() {

        if (!likelihoodKnown) {
            logLikelihood = calculateLogLikelihood();
            likelihoodKnown = true;
        }
        return logLikelihood;
    }

    private double calculateLogLikelihood() {

        double logLike = 0.0;

        for (int i = 0; i < tree.getInternalNodeCount(); ++i) {
            NodeRef node = tree.getNode(i);
            if (node != tree.getRoot()) {
                if (updateNeededForNode(tree, node)) {
                    double branchLength = tree.getBranchLength(node);
                    double latentProportion = getLatentProportion(tree, node);

                    assert(latentProportion < 1.0);

                    double reward = branchLength * latentProportion;
                    double density = getBranchRewardDensity(reward, branchLength);
                    branchLikelihoods[node.getNumber()] = Math.log(density);
                }
                logLike += branchLikelihoods[node.getNumber()];
                // TODO More importantly, MH proposals on [0,1] may be missing a Jacobian for which we should adjust.
                // TODO This is easy to test and we should do it when sampling appears to work.
            }
        }

        clearUpdateAllBranches();
        clearAllCategories();

        return logLike;
    }

    private boolean updateNeededForNode(Tree tree, NodeRef node) {
        if (USE_CACHING) {
            return (updateCategory != null && updateCategory[branchCategoryProvider.getBranchCategory(tree, node)]) || updateBranch[node.getNumber()];
        } else {
            return true;
        }
    }

    public double getBranchRewardDensity(double reward, double branchLength) {
        if (series == null) {
            series = createSeries();
        }

        int state = 0 * 2 + 0; // just start = end = 0 entry
        // Reward is [0,1], and we want to track time in latent state (= 1).
        // Therefore all nodes are in state 0
//        double joint = series.computePdf(reward, branchLength)[state];
        double joint = series.computePdf(reward, branchLength, 0, 0);
        double marg = series.computeConditionalProbability(branchLength, 0, 0);
        // TODO Overhead in creating double[] could be saved by changing signature to computePdf

        return joint / marg; // conditional on ending state.
    }

    @Override
    public void makeDirty() {
        likelihoodKnown = false;
        series = null;
        setUpdateAllBranches();
    }

    @Override
    public String getTraitName() {
        return BranchRateModel.RATE;
    }

    @Override
    public Intent getIntent() {
        return Intent.BRANCH;
    }

    @Override
    public TreeTrait getTreeTrait(final String key) {
        if (key.equals(BranchRateModel.RATE)) {
            return this;
        } else if (latentStateProportions != null && key.equals(latentStateProportions.getTraitName())) {
            return latentStateProportions;
        } else if (branchCategoryProvider != null && key.equals(branchCategoryProvider.getTraitName())) {
            return branchCategoryProvider;
        } else {
            throw new IllegalArgumentException("Unrecognised Tree Trait key, " + key);
        }
    }

    @Override
    public TreeTrait[] getTreeTraits() {
        return new TreeTrait[]{this, latentStateProportions, branchCategoryProvider};
    }

    @Override
    public Class getTraitClass() {
        return Double.class;
    }

    @Override
    public boolean getLoggable() {
        return true;
    }

    @Override
    public Double getTrait(final Tree tree, final NodeRef node) {
        return getBranchRate(tree, node);
    }

    @Override
    public String getTraitString(final Tree tree, final NodeRef node) {
        return Double.toString(getBranchRate(tree, node));
    }

    static class Mode {
        double pdf;
        double reward;

        Mode(double pdf, double reward) {
            this.pdf = pdf;
            this.reward = reward;
        }
    }

    static Mode findMode(List<Double> values, List<Double> rewards) {
        Mode find = new Mode(values.get(0), rewards.get(0));
        for (int i = 1; i < values.size(); ++i) {
            if (values.get(i) > find.pdf) {
                find.pdf = values.get(i);
                find.reward = rewards.get(i);
            }
        }
        return find;
    }

    static double calculateExpectation(List<Double> pdfs, List<Double> rewards) {
        double weight = 0.0;
        double wsum = 0.0;
        for (int i = 0; i < pdfs.size(); ++i) {
            weight += pdfs.get(i);
            wsum += rewards.get(i) * pdfs.get(i);
        }
        double mean = wsum / weight;
//        System.err.println(wsum);
//        System.err.println(weight);
//        System.err.println(mean);
//
//        System.exit(-1);
        return wsum;
    }

    public static void main(String[] args) {

        Parameter rate = new Parameter.Default(2.0);
        Parameter prop = new Parameter.Default(0.5);


        SericolaLatentStateBranchRateModel model = new SericolaLatentStateBranchRateModel(rate, prop);


        for (double branchLength = 0.1; branchLength <= 10.0; branchLength += 0.1) {

            List<Double> pdfs = new ArrayList<Double>();
            List<Double> rewards = new ArrayList<Double>();

            for (double reward = 0; reward <= branchLength; reward += 0.01 * branchLength) {
                double value = model.getBranchRewardDensity(reward, branchLength);
//                System.out.println(reward + "," + model.getBranchRewardDensity(reward, branchLength));
                System.out.println();
                rewards.add(reward);
                pdfs.add(value);
            }
            Mode mode = findMode(pdfs, rewards);

//        System.out.println();
            System.out.println(branchLength // ", " + mode.reward //+ " " + mode.pdf
                    + " " + (mode.reward / branchLength)
                    + " " + (calculateExpectation(pdfs, rewards) / branchLength)
            );
            //System.out.println(model.getSeries());
        }

    }

    public MarkovReward getSeries() {
        if (series == null) {
            series = createSeries();
        }
        return series;
    }
}
