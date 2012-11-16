/*
 * EpochBranchModel.java
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

package dr.app.beagle.evomodel.branchmodel;

import java.util.ArrayList;
import java.util.List;

import dr.app.beagle.evomodel.substmodel.SubstitutionModel;
import dr.evolution.tree.NodeRef;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.util.Author;
import dr.util.Citable;
import dr.util.Citation;

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
        int lastTransitionTime = nModels - 2;

        double[] transitionTimes = epochTimes.getParameterValues();

        double parentHeight = tree.getNodeHeight(tree.getParent(node));
        double nodeHeight = tree.getNodeHeight(node);
        double branchLength = tree.getBranchLength(node);

        List<Double> weightList = new ArrayList<Double>();
        List<Integer> orderList = new ArrayList<Integer>();

        if (parentHeight <= transitionTimes[0]) {

            weightList.add( branchLength );
            orderList.add(0);

        } else {

            // first case: 0-th transition time
            if (nodeHeight < transitionTimes[0] && transitionTimes[0] <= parentHeight) {

                weightList.add( transitionTimes[0] - nodeHeight );
                orderList.add(0);

            } else {
                // do nothing
            }// END: 0-th model check

            // second case: i to i+1 transition times
            for (int i = 1; i <= lastTransitionTime; i++) {

                if (nodeHeight < transitionTimes[i]) {

                    if (parentHeight <= transitionTimes[i] && transitionTimes[i - 1] < nodeHeight) {

                        weightList.add( branchLength );
                        orderList.add(i);

                    } else {

                        double startTime = Math.max(nodeHeight, transitionTimes[i - 1]);
                        double endTime = Math.min(parentHeight, transitionTimes[i]);

                        if (endTime >= startTime) {

                            weightList.add( endTime - startTime );
                            orderList.add(i);

                        }// END: negative weights check

                    }// END: full branch in middle epoch check

                }// END: i-th model check

            }// END: i loop

            // third case: last transition time
            if (parentHeight >= transitionTimes[lastTransitionTime] && transitionTimes[lastTransitionTime] > nodeHeight) {

                weightList.add( parentHeight - transitionTimes[lastTransitionTime] );
                orderList.add(nModels - 1);

            } else if (nodeHeight > transitionTimes[lastTransitionTime]) {

                weightList.add( branchLength );
                orderList.add(nModels - 1);

            } else {
                // nothing to add
            }// END: last transition time check

        }// END: if branch below first transition time bail out

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


    /**
     * @return a list of citations associated with this object
     */
    public List<Citation> getCitations() {
        List<Citation> citations = new ArrayList<Citation>();
        citations.add(new Citation(new Author[]{new Author("F", "Bielejec"),
                new Author("P", "Lemey"), new Author("G", "Baele"), new Author("A", "Rambaut"),
                new Author("MA", "Suchard")}, Citation.Status.IN_PREPARATION));
        return citations;
    }// END: getCitations

    private final TreeModel tree;
    private final List<SubstitutionModel> substitutionModels;
    private final Parameter epochTimes;
}// END: class
