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

import dr.app.beagle.evomodel.substmodel.FrequencyModel;
import dr.app.beagle.evomodel.substmodel.GeneralSubstitutionModel;
import dr.app.beagle.evomodel.substmodel.SubstitutionModel;
import dr.app.beagle.evomodel.substmodel.UniformizedSubstitutionModel;
import dr.evolution.datatype.TwoStates;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.tree.TreeParameterModel;
import dr.inference.markovjumps.MarkovJumpsType;
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
                                              Parameter latentStateProportionParameter) {
        this(LATENT_STATE_BRANCH_RATE_MODEL, treeModel, nonLatentRateModel, latentTransitionRateParameter, latentStateProportionParameter);
    }

    public SericolaLatentStateBranchRateModel(String name,
                                              TreeModel treeModel,
                                              BranchRateModel nonLatentRateModel,
                                              Parameter latentTransitionRateParameter,
                                              Parameter latentStateProportionParameter) {
        super(name);

        this.tree = treeModel;
        addModel(tree);

        this.nonLatentRateModel = nonLatentRateModel;
        addModel(nonLatentRateModel);

        this.latentTransitionRateParameter = latentTransitionRateParameter;
        addVariable(latentTransitionRateParameter);

        this.latentStateProportionParameter = latentStateProportionParameter;
        addVariable(latentTransitionRateParameter);

        this.latentStateProportions = new TreeParameterModel(tree, latentStateProportionParameter, false);
    }

    /**
     * Empty constructor for use by the main
     */
    public SericolaLatentStateBranchRateModel() {
        super(LATENT_STATE_BRANCH_RATE_MODEL);
        tree = null;
        nonLatentRateModel = null;
        latentTransitionRateParameter = null;
        latentStateProportionParameter = null;
        latentStateProportions = null;

    }

    private double[] createLatentInfinitesimalMatrix() {
        final double rate = latentTransitionRateParameter.getParameterValue(0);
        final double prop = latentStateProportionParameter.getParameterValue(0);

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
        return new SericolaSeriesMarkovReward(createLatentInfinitesimalMatrix(), createReward(), 2);
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

        } else if (model == nonLatentRateModel) {

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
        if (variable == latentStateProportionParameter || variable == latentTransitionRateParameter) {
            // series computations have changed
            series = null;
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

        if (series == null) {
            series = createSeries();
        }

        double logLike = 0.0;
        for (int i = 0; i < tree.getInternalNodeCount(); ++i) {
            NodeRef node = tree.getNode(i);
            if (node != tree.getRoot()) {
                double branchLength = tree.getBranchLength(node);
                double latentProportion = latentStateProportions.getNodeValue(tree, node);
                double reward = branchLength * latentProportion;

                logLike += Math.log(series.computePdf(reward, branchLength)[1 * 2 + 1]); // just start = end = 1
                // TODO These quantities are off by the log(conditional transition probability) alone each branch,
                // TODO but cprobs should be constant with Metropolis-Hastings proposal on proportions alone
            }
        }

        return logLike;
    }

    @Override
    public void makeDirty() {

    }

    @Override
    public String getTraitName() {
        return null;
    }

    @Override
    public Intent getIntent() {
        return null;
    }

    @Override
    public Class getTraitClass() {
        return null;
    }

    @Override
    public Double getTrait(Tree tree, NodeRef node) {
        return null;
    }

    @Override
    public String getTraitString(Tree tree, NodeRef node) {
        return null;
    }

    @Override
    public boolean getLoggable() {
        return false;
    }

    @Override
    public TreeTrait[] getTreeTraits() {
        return new TreeTrait[0];
    }

    @Override
    public TreeTrait getTreeTrait(String key) {
        return null;
    }

    public static void main(String[] args) {
        FrequencyModel frequencyModel = new FrequencyModel(TwoStates.INSTANCE, new double[]{0.5, 0.5});
        Parameter rateParameter = new Parameter.Default(new double[]{0.001, 0.001});

        SubstitutionModel binaryModel = new GeneralSubstitutionModel("binary", TwoStates.INSTANCE, frequencyModel, rateParameter, 0);

        // Do this once, for all branches
        UniformizedSubstitutionModel uSM = new UniformizedSubstitutionModel(binaryModel, MarkovJumpsType.REWARDS);
        uSM.setSaveCompleteHistory(true);
        double[] rewardRegister = new double[]{0.0, 1.0};
        uSM.setRegistration(rewardRegister);

        double delta = 1;
        double[] values = new double[101];
        int count = 100000;

        double length = 0.0;
        for (int j = 0; j < 100; j++) {
            double reward = uSM.computeCondStatMarkovJumps(0, 0, length);
            double proportionTime = reward / length;
            System.out.println(length + "\t" + proportionTime);
            length += delta;
        }

        // average over many realizations...
//        for (int i = 0; i < count; i++) {
//            double length = 0.0;
//            for (int j = 0; j < 100; j++) {
//                double reward = uSM.computeCondStatMarkovJumps(0, 0, length);
//                double proportionTime = reward / length;
//                values[j] += proportionTime;
//                length += delta;
//            }
//        }
//
//        double length = 0.0;
//        for (int j = 0; j < 100; j++) {
//            System.out.println(length + "\t" + values[j] / count);
//            length += delta;
//        }
    }
}
