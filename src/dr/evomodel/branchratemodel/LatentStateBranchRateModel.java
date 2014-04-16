/*
 * LatentStateBranchRateModel.java
 *
 * Copyright (C) 2002-2014 Alexei Drummond, Andrew Rambaut & Marc Suchard
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
import dr.inference.model.AbstractModelLikelihood;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;

/**
 * LatentStateBranchRateModel
 *
 * @author Andrew Rambaut
 * @author Marc Suchard
 * @version $Id$
 *
 * $HeadURL$
 *
 * $LastChangedBy$
 * $LastChangedDate$
 * $LastChangedRevision$
 */
public class LatentStateBranchRateModel extends AbstractModelLikelihood implements BranchRateModel {

    public static final String LATENT_STATE_BRANCH_RATE_MODEL = "latentStateBranchRateModel";

    private final TreeModel tree;
    private final BranchRateModel nonLatentRateModel;
    private final Parameter latentTransitionRateParameter;
    private final Parameter latentStateProportionParameter;
    private final TreeParameterModel latentStateProportions;

    private UniformizedSubstitutionModel uSM = null;

    public LatentStateBranchRateModel(TreeModel treeModel,
                                      BranchRateModel nonLatentRateModel,
                                      Parameter latentTransitionRateParameter,
                                      Parameter latentStateProportionParameter) {
        this(LATENT_STATE_BRANCH_RATE_MODEL, treeModel, nonLatentRateModel, latentTransitionRateParameter, latentStateProportionParameter);
    }

    public LatentStateBranchRateModel(String name,
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

        if (latentStateProportionParameter != null) {
            this.latentStateProportionParameter = latentStateProportionParameter;
            addVariable(latentTransitionRateParameter);
        } else {
            this.latentStateProportionParameter = new Parameter.Default(0.5);
        }
        this.latentStateProportions = new TreeParameterModel(tree, latentStateProportionParameter, false);
    }

    /**
     * Empty constructor for use by the main
     */
    public LatentStateBranchRateModel() {
        super(LATENT_STATE_BRANCH_RATE_MODEL);
        tree = null;
        nonLatentRateModel = null;
        latentTransitionRateParameter = null;
        latentStateProportionParameter = null;
        latentStateProportions = null;

    }

    @Override
    public double getBranchRate(Tree tree, NodeRef node) {
        if (uSM == null) {
            createRewards();
        }

        double length = tree.getBranchLength(node);

        double nonLatentRate = nonLatentRateModel.getBranchRate(tree, node);

        double latentProportion;
        if (latentStateProportionParameter != null) {
            latentProportion = latentStateProportions.getNodeValue(tree, node);
        } else {
            latentProportion = calculateLatentProportion(length);
            latentStateProportions.setNodeValue(tree, node, latentProportion);
        }

        return calculateBranchRate(nonLatentRate, latentProportion);
    }

    public double getLatentProportion(double length) {
        if (uSM == null) {
            createRewards();
        }

        return calculateLatentProportion(length);
    }

    private double calculateLatentProportion(double length) {
        // Do this each time you need a new branch rate
        double reward = uSM.computeCondStatMarkovJumps(0, 0, length);
        double proportionTime = reward / length;
        return proportionTime;
    }

    private double calculateBranchRate(double nonLatentRate, double latentProportion) {
        return nonLatentRate * (1.0 - latentProportion);
    }

    private void createRewards() {
        FrequencyModel frequencyModel = new FrequencyModel(TwoStates.INSTANCE, new double[] { 0.5, 0.5 });
        Parameter rateParameter = new Parameter.Default(new double[] { 1.0, 1.0 });

        SubstitutionModel binaryModel = new GeneralSubstitutionModel("binary", TwoStates.INSTANCE, frequencyModel, rateParameter, 0);

        // Do this once, for all branches
        uSM = new UniformizedSubstitutionModel(binaryModel, MarkovJumpsType.REWARDS);
        uSM.setSaveCompleteHistory(true);
        double[] rewardRegister = new double[]{0.0, 1.0};
        uSM.setRegistration(rewardRegister);
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        if (model == tree) {

        } else if (model == nonLatentRateModel) {

        }

    }

    @Override
    protected void storeState() {

    }

    @Override
    protected void restoreState() {

    }

    @Override
    protected void acceptState() {

    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {

    }

    @Override
    public Model getModel() {
        return this;
    }

    @Override
    public double getLogLikelihood() {
        throw new UnsupportedOperationException("Not implemented yet");
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
        LatentStateBranchRateModel lsbrm = new LatentStateBranchRateModel();
        double delta = 0.01;
        double[] values = new double[101];
        int count = 100000;

        for (int i = 0; i < count; i++) {
            double length = 0.01;
            for (int j = 0; j < 100; j++) {
                values[j] += lsbrm.getLatentProportion(length);
            }
        }
        double length = 0.01;
        for (int j = 0; j < 100; j++) {
            System.out.println(length + "\t" + values[j] / count);
            length += delta;
        }
    }
}
