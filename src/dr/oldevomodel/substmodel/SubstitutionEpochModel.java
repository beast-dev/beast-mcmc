/*
 * SubstitutionEpochModel.java
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

package dr.oldevomodel.substmodel;

import dr.evolution.datatype.DataType;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.math.matrixAlgebra.Vector;

import java.util.List;

/**
 * @author Marc A. Suchard
 */

public class SubstitutionEpochModel extends AbstractSubstitutionModel {

    public static final boolean DEBUG = false;

    public SubstitutionEpochModel(String name,
                           List<SubstitutionModel> modelList,
                           Parameter transitionTimes,
                           DataType dataType, FrequencyModel freqModel) {

        super(name, dataType, freqModel);

        this.modelList = modelList;
        this.transitionTimesParameter = transitionTimes;
        this.transitionTimes = transitionTimesParameter.getParameterValues();

        addVariable(transitionTimes);

        for (SubstitutionModel model : modelList)
            addModel(model);

        numberModels = modelList.size();
        weight = new double[numberModels];
        stateCount = dataType.getStateCount();
        stepMatrix = new double[stateCount * stateCount];
        productMatrix = new double[stateCount * stateCount];
        resultMatrix = new double[stateCount * stateCount];
    }


    protected void frequenciesChanged() {
    }

    protected void handleModelChangedEvent(Model model, Object object, int index) {
        if (model == freqModel)
            frequenciesChanged();
        else // This is an epoch model and I need to pass the info on
            fireModelChanged(object, index);
    }

    protected void ratesChanged() {
    }

    protected void setupRelativeRates() {
    }

    public void getTransitionProbabilities(double startTime, double endTime, double distance, double[] matrix) {
        int matrixCount = 0;
        
//        System.out.println("startTime " + startTime ); 
//        System.out.println("endTime " + endTime );  
//            EpochBranchSubstitutionModel.printMatrix(resultMatrix);
        
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
                    
//                    System.out.println("first " + weight[m] * (endTime - startTime) + " " + model.getVariable(0).getValue(0) + model.getId());
//                    EpochBranchSubstitutionModel.printMatrix(resultMatrix);
                    
                } else {
                    model.getTransitionProbabilities(distance * weight[m], stepMatrix);
                    // Sum over unobserved state
                    int index = 0;
                    
//                    System.out.println("startTime " + startTime ); 
//                    System.out.println("endTime " + endTime );  
//                    System.out.println("second " + weight[m] * (endTime - startTime) + " " + model.getVariable(0).getValue(0) + model.getId());
//                        EpochBranchSubstitutionModel.printMatrix(stepMatrix);
                    
                    for (int i = 0; i < stateCount; i++) {
                        for (int j = 0; j < stateCount; j++) {
                            productMatrix[index] = 0;
                            for (int k = 0; k < stateCount; k++) {
                                productMatrix[index] += resultMatrix[i * stateCount + k] * stepMatrix[k * stateCount + j];
                            }
                            index++;
                        }
                    }
                 
//            		EpochBranchSubstitutionModel.printMatrix(productMatrix);
                    
                    // Swap pointers
                    double[] tmpMatrix = resultMatrix;
                    resultMatrix = productMatrix;
                    productMatrix = tmpMatrix;
                }
            }
        }
        if (!oneMatrix)
            System.arraycopy(productMatrix, 0, matrix, 0, stateCount * stateCount);
        
//		System.out.println("C:");
//		EpochBranchSubstitutionModel.printMatrix(resultMatrix);
//        System.exit(-1);

    }

    private int getEpochWeights(double startTime, double endTime, double[] weights) {

        int matrixCount = 0;
        final double lengthTime = endTime - startTime;
        final int lastTime = numberModels - 2;

        // model 0, 1, 2, ..., K-2, K-1
        // times   0, 1,  ...,   K-2,
        // where K = numberModels

//        System.out.println(lengthTime);
        
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

//        System.out.println(endTime-startTime);
//        System.out.println(matrixCount);
//        EpochBranchSubstitutionModel.printArray(weights);
        	
        return matrixCount;
    }

    public void getTransitionProbabilities(double distance, double[] matrix) {
        throw new RuntimeException("Should not get here in a substitution epoch model.");
    }

    protected final void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        super.handleVariableChangedEvent(variable, index, type);

        if (variable == transitionTimesParameter) {
            transitionTimes = transitionTimesParameter.getParameterValues();
            fireModelChanged(variable, index);
        }
    }

    protected void storeState() {
    }

    protected void restoreState() {
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
