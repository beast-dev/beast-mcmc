/*
 * EpochBranchModel.java
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

package dr.evomodel.branchmodel;

import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evolution.tree.NodeRef;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.util.Author;
import dr.util.Citable;
import dr.util.Citation;

import java.util.*;

/**
 * @author Filip Bielejec
 * @author Andrew Rambaut
 * @author Marc A. Suchard
 * @version $Id$
 */
@SuppressWarnings("serial")
public class EpochBranchModel extends AbstractModel implements BranchModel, Citable {

    public static final String EPOCH_BRANCH_MODEL = "EpochBranchModel";


    public EpochBranchModel(TreeModel tree,
                            List<SubstitutionModel> substitutionModels,
                            Parameter epochTimes) {

        super(EPOCH_BRANCH_MODEL);

        this.substitutionModels = substitutionModels;
        if (substitutionModels == null || substitutionModels.size() == 0) {
            throw new IllegalArgumentException("EpochBranchModel must be provided with at least one substitution model");
        }

        this.epochTimes = epochTimes;
        this.tree = tree;

        for (SubstitutionModel model : substitutionModels) {
            addModel(model);
        }

        addModel(tree);
        addVariable(epochTimes);
    }// END: Constructor

    @Override
    public Mapping getBranchModelMapping(NodeRef node) {

        int nModels = substitutionModels.size();
        int epochCount = nModels - 1;

        double[] transitionTimes = epochTimes.getParameterValues();

        double parentHeight = tree.getNodeHeight(tree.getParent(node));
        double nodeHeight = tree.getNodeHeight(node);

        List<Double> weightList = new ArrayList<Double>();
        List<Integer> orderList = new ArrayList<Integer>();

        // find the epoch that the node height is in...
        int epoch = 0;
        while (epoch < epochCount && nodeHeight >= transitionTimes[epoch]) {
            epoch ++;
        }

        double currentHeight = nodeHeight;

        // find the epoch that the parent height is in...
        while (epoch < epochCount && parentHeight >= transitionTimes[epoch]) {
            weightList.add( transitionTimes[epoch] - currentHeight );
            orderList.add(epoch);

            currentHeight = transitionTimes[epoch];

            epoch ++;
        }

        weightList.add( parentHeight - currentHeight );
        orderList.add(epoch);


        if (orderList.size() == 0) {
            throw new RuntimeException("EpochBranchModel failed to give a valid mapping");
        }

        final int[] order = new int[orderList.size()];
        final double[] weights = new double[weightList.size()];
        for (int i = 0; i < orderList.size(); i++) {
            order[i] = orderList.get(i);
            weights[i] = weightList.get(i);
        }

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
    }// END: getBranchModelMapping

    @Override
    public boolean requiresMatrixConvolution() {
        return true;
    }

    @Override
    public List<SubstitutionModel> getSubstitutionModels() {
        return substitutionModels;
    }

    @Override
    public SubstitutionModel getRootSubstitutionModel() {
        return substitutionModels.get(substitutionModels.size() - 1);
    }

    public FrequencyModel getRootFrequencyModel() {
        return getRootSubstitutionModel().getFrequencyModel();
    }

    protected void handleModelChangedEvent(Model model, Object object, int index) {
        fireModelChanged();
    }// END: handleModelChangedEvent

    @SuppressWarnings("rawtypes")
    protected void handleVariableChangedEvent(Variable variable, int index,
                                              Parameter.ChangeType type) {
    }// END: handleVariableChangedEvent

    protected void storeState() {
    }// END: storeState

    protected void restoreState() {
    }// END: restoreState

    protected void acceptState() {
    }// END: acceptState

    @Override
    public Citation.Category getCategory() {
        return Citation.Category.SUBSTITUTION_MODELS;
    }

    @Override
    public String getDescription() {
        return "Epoch Branch model";
    }

    public List<Citation> getCitations() {
        return Arrays.asList(
                new Citation(new Author[]{new Author("F", "Bielejec"),
                        new Author("P", "Lemey"), new Author("G", "Baele"), new Author("A", "Rambaut"),
                        new Author("MA", "Suchard")}, Citation.Status.IN_PREPARATION));
    }// END: getCitations

    private final TreeModel tree;
    private final List<SubstitutionModel> substitutionModels;
    private final Parameter epochTimes;
}// END: class
