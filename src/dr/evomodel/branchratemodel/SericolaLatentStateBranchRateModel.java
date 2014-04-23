/*
 * SericolaLatentStateBranchRateModel.java
 *
 * Copyright (c) 2002-2014 Alexei Drummond, Andrew Rambaut and Marc Suchard
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
import dr.inference.markovjumps.SericolaSeriesMarkovReward;
import dr.inference.model.AbstractModelLikelihood;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;

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

    private final TreeModel tree;
    private final BranchRateModel nonLatentRateModel;
    private final Parameter latentTransitionRateParameter;
    private final Parameter latentTransitionFrequencyParameter;
    private final Parameter latentStateProportionParameter;
    private final TreeParameterModel latentStateProportions;

    private SericolaSeriesMarkovReward series;
    private SericolaSeriesMarkovReward storedSeries;
    private boolean likelihoodKnown = false;
    private boolean storedLikelihoodKnown;
    private double logLikelihood;
    private double storedLogLikelihood;

//    private UniformizedSubstitutionModel uSM = null;

    public SericolaLatentStateBranchRateModel(TreeModel treeModel,
                                              BranchRateModel nonLatentRateModel,
                                              Parameter latentTransitionRateParameter,
                                              Parameter latentTransitionFrequencyParameter,
                                              Parameter latentStateProportionParameter) {
        this(LATENT_STATE_BRANCH_RATE_MODEL, treeModel, nonLatentRateModel, latentTransitionRateParameter,
                latentTransitionFrequencyParameter, latentStateProportionParameter);
    }

    public SericolaLatentStateBranchRateModel(String name,
                                              TreeModel treeModel,
                                              BranchRateModel nonLatentRateModel,
                                              Parameter latentTransitionRateParameter,
                                              Parameter latentTransitionFrequencyParameter,
                                              Parameter latentStateProportionParameter) {
        super(name);

        this.tree = treeModel;
        addModel(tree);

        this.nonLatentRateModel = nonLatentRateModel;
        addModel(nonLatentRateModel);

        this.latentTransitionRateParameter = latentTransitionRateParameter;
        addVariable(latentTransitionRateParameter);

        this.latentTransitionFrequencyParameter = latentTransitionFrequencyParameter;
        addVariable(latentTransitionFrequencyParameter);

        this.latentStateProportionParameter = latentStateProportionParameter;
        addVariable(latentStateProportionParameter);   // TODO This may not be necessary

        this.latentStateProportions = new TreeParameterModel(tree, latentStateProportionParameter, false, Intent.BRANCH);
        addModel(latentStateProportions);
    }

    /**
     * Empty constructor for use by the main
     */
    public SericolaLatentStateBranchRateModel() {
        super(LATENT_STATE_BRANCH_RATE_MODEL);
        tree = null;
        nonLatentRateModel = null;
        latentTransitionRateParameter = null;
        latentTransitionFrequencyParameter = null;
        latentStateProportionParameter = null;
        latentStateProportions = null;

    }

    public SericolaLatentStateBranchRateModel(Parameter rate, Parameter prop) {
        super(LATENT_STATE_BRANCH_RATE_MODEL);
        tree = null;
        nonLatentRateModel = null;
        latentTransitionRateParameter = rate;
        latentTransitionFrequencyParameter = prop;
        latentStateProportionParameter = null;
        latentStateProportions = null;
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

    private SericolaSeriesMarkovReward createSeries() {
        SericolaSeriesMarkovReward series = new SericolaSeriesMarkovReward(createLatentInfinitesimalMatrix(),
                createReward(), 2);
        return series;
    }

    @Override
    public double getBranchRate(Tree tree, NodeRef node) {

        double nonLatentRate = nonLatentRateModel.getBranchRate(tree, node);
        double latentProportion = latentStateProportions.getNodeValue(tree, node);

        return calculateBranchRate(nonLatentRate, latentProportion);
    }

    private double calculateBranchRate(double nonLatentRate, double latentProportion) {
        return nonLatentRate * (1.0 - latentProportion);
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        if (model == tree) {
            likelihoodKnown = false; // node heights change elasped times on branches, TODO could cache
        } else if (model == nonLatentRateModel) {
            // rates will change but the latent proportions haven't so the density is unchanged
            fireModelChanged();
        } else if (model == latentStateProportions) {
            fireModelChanged();
            likelihoodKnown = false; // argument of density has changed
        }
    }

    @Override
    protected void storeState() {
        storedSeries = series;
        storedLogLikelihood = logLikelihood;
        storedLikelihoodKnown = likelihoodKnown;
    }

    @Override
    protected void restoreState() {
        series = storedSeries;
        logLikelihood = storedLogLikelihood;
        likelihoodKnown = storedLikelihoodKnown;
    }

    @Override
    protected void acceptState() {

    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        if (variable == latentTransitionFrequencyParameter || variable == latentTransitionRateParameter) {
            // series computations have changed
            series = null;
            likelihoodKnown = false;
        } else if (variable == latentStateProportionParameter) {
            likelihoodKnown = false;
        }
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
                double branchLength = tree.getBranchLength(node);
                double latentProportion = latentStateProportions.getNodeValue(tree, node);
                double reward = branchLength * latentProportion;
                double density = getBranchRewardDensity(reward, branchLength);
                logLike += Math.log(density);
                // TODO More importantly, MH proposals on [0,1] may be missing a Jacobian for which we should adjust.
                // TODO This is easy to test and we should do it when sampling appears to work.
            }
        }

        return logLike;
    }

    public double getBranchRewardDensity(double reward, double branchLength) {
        if (series == null) {
            series = createSeries();
        }
        int state = 0 * 2 + 0; // just start = end = 0 entry
        // Reward is [0,1], and we want to track time in latent state (= 1).
        // Therefore all nodes are in state 0
        double joint = series.computePdf(reward, branchLength)[state];
        double marg = series.computeConditionalProbability(branchLength, 0, 0);
        // TODO Overhead in creating double[] could be saved by changing signature to computePdf

        return joint / marg; // conditional on ending state.
    }

    @Override
    public void makeDirty() {
        likelihoodKnown = false;
        series = null;
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
        } else if (key.equals(latentStateProportions.getTraitName())) {
            return latentStateProportions;
        } else {
            throw new IllegalArgumentException("Unrecognised Tree Trait key, " + key);
        }
    }

    @Override
    public TreeTrait[] getTreeTraits() {
        return new TreeTrait[]{this, latentStateProportions};
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

    public static void main(String[] args) {

        Parameter rate = new Parameter.Default(4.4);
        Parameter prop = new Parameter.Default(0.25);

        SericolaLatentStateBranchRateModel model = new SericolaLatentStateBranchRateModel(rate, prop);

        double branchLength = 2.0;
        for (double reward = 0; reward < branchLength; reward += 0.01) {
            System.out.println(reward + ",\t" + model.getBranchRewardDensity(reward, branchLength) + ",");
        }

        System.out.println();
        System.out.println(model.getSeries());

    }

    public SericolaSeriesMarkovReward getSeries() {
        if (series == null) {
            series = createSeries();
        }
        return series;
    }
}
