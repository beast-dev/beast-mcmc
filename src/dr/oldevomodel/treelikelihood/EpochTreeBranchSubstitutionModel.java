/*
 * EpochTreeBranchSubstitutionModel.java
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

package dr.oldevomodel.treelikelihood;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.oldevomodel.sitemodel.SiteModel;
import dr.oldevomodel.substmodel.SubstitutionModel;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.math.matrixAlgebra.Vector;

import java.util.List;

/**
 * @author Marc A. Suchard
 */
@Deprecated // Switching to BEAGLE
public class EpochTreeBranchSubstitutionModel extends TreeBranchSubstitutionModel {

    public static final boolean DEBUG = false;

    public EpochTreeBranchSubstitutionModel(String name,
                                            SiteModel siteModel, List<SubstitutionModel> substModelList, BranchRateModel branchModel,
                                            Parameter transitionTimes) {
        super(name, siteModel, null, branchModel);
        this.modelList = substModelList;

        this.transitionTimesParameter = transitionTimes;
        this.transitionTimes = transitionTimesParameter.getParameterValues();

        addVariable(transitionTimes);

        for (SubstitutionModel model : modelList)
            addModel(model);

        numberModels = modelList.size();
        weight = new double[numberModels];
        stateCount = modelList.get(0).getDataType().getStateCount();
        stepMatrix = new double[stateCount * stateCount];
        productMatrix = new double[stateCount * stateCount];
        resultMatrix = new double[stateCount * stateCount];
    }

    public void getTransitionProbabilities(Tree tree, NodeRef node, int rateCategory, double[] matrix) {

        NodeRef parent = tree.getParent(node);

        final double branchRate = branchModel.getBranchRate(tree, node);

        // Get the operational time of the branch
        final double startTime = tree.getNodeHeight(parent);
        final double endTime = tree.getNodeHeight(node);
        final double branchTime = branchRate * (startTime - endTime);

        if (branchTime < 0.0) {
            throw new RuntimeException("Negative branch length: " + branchTime);
        }

        double distance = siteModel.getRateForCategory(rateCategory) * branchTime;

        int matrixCount = 0;
        boolean oneMatrix = (getEpochWeights(startTime, endTime, weight) == 1);
        for (int m = 0; m < numberModels; m++) {
            if (weight[m] > 0) {
                SubstitutionModel model = modelList.get(m);
                if (matrixCount == 0) {
                    if (oneMatrix) {
                        model.getTransitionProbabilities(distance, matrix);
                        break;
                    } else
                        model.getTransitionProbabilities(distance * weight[m], resultMatrix);
                    matrixCount++;
                } else {
                    model.getTransitionProbabilities(distance * weight[m], stepMatrix);
                    // Sum over unobserved state
                    int index = 0;
                    for (int i = 0; i < stateCount; i++) {
                        for (int j = 0; j < stateCount; j++) {
                            productMatrix[index] = 0;
                            for (int k = 0; k < stateCount; k++) {
                                productMatrix[index] += resultMatrix[i * stateCount + k] * stepMatrix[k * stateCount + j];
                            }
                            index++;
                        }
                    }
                    // Swap pointers
                    double[] tmpMatrix = resultMatrix;
                    resultMatrix = productMatrix;
                    productMatrix = tmpMatrix;
                }
            }
        }
        if (!oneMatrix)
            System.arraycopy(resultMatrix, 0, matrix, 0, stateCount * stateCount);
    }

    private int getEpochWeights(double startTime, double endTime, double[] weights) {

        int matrixCount = 0;
        final double lengthTime = endTime - startTime;
        final int lastTime = numberModels - 2;

        // model 0, 1, 2, ..., K-2, K-1
        // times   0, 1,  ...,   K-2,
        // where K = numberModels

        // First epoch: 0 -> transitionTimes[0];
        if (startTime <= transitionTimes[0]) {
            if (endTime <= transitionTimes[0])
                weights[0] = 1;
            else
                weights[0] = (transitionTimes[0] - startTime) / lengthTime;
            matrixCount++;
        } else
            weights[0] = 0;

        // Middle epoches:
        for (int i = 1; i <= lastTime; i++) {
            if (startTime <= transitionTimes[i]) {
                double start = Math.max(startTime, transitionTimes[i - 1]);
                double end = Math.min(endTime, transitionTimes[i]);
                weights[i] = (end - start) / lengthTime;
                matrixCount++;
            } else
                weights[i] = 0;
        }

        // Last epoch: transitionTimes[K-2] -> Infinity
        if (lastTime >= 0) {
            if (endTime > transitionTimes[lastTime]) {
                double start = Math.max(startTime, transitionTimes[lastTime]);
                weights[lastTime + 1] = (endTime - start) / lengthTime;
                matrixCount++;
            } else
                weights[lastTime + 1] = 0;
        }

        if (DEBUG) {
            double totalWeight = 0;
            for (int i = 0; i < numberModels; i++)
                totalWeight += weights[i];
            System.err.println("Start: " + startTime + " End: " + endTime + " Count: " + matrixCount + " Weight: " + totalWeight + " - " + new Vector(weights));
            if (totalWeight > 1.001) System.exit(-1);
            if (totalWeight < 0.999) System.exit(-1);
        }

        return matrixCount;
    }

    protected final void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        super.handleVariableChangedEvent(variable, index, type);

        if (variable == transitionTimesParameter) {
            transitionTimes = transitionTimesParameter.getParameterValues();
            fireModelChanged(variable, index);
        }
    }

    private List<SubstitutionModel> modelList;
    private Parameter transitionTimesParameter;
    private double[] transitionTimes;
    private double[] weight;
    private double[] stepMatrix;
    private double[] productMatrix;
    private double[] resultMatrix;
    private int numberModels;
    private int stateCount;

}

