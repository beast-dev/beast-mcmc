/*
 * MarkovJumpsSubstitutionModel.java
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

package dr.evomodel.substmodel;

import dr.evolution.datatype.DataType;
import dr.inference.markovjumps.MarkovJumpsCore;
import dr.inference.markovjumps.MarkovJumpsType;
import dr.inference.markovjumps.StateHistory;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;

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

public class MarkovJumpsSubstitutionModel extends AbstractModel {

    public MarkovJumpsSubstitutionModel(SubstitutionModel substModel) {
        this(substModel, MarkovJumpsType.COUNTS);
    }

    public MarkovJumpsSubstitutionModel(SubstitutionModel substModel, MarkovJumpsType type) {
        super(substModel.getModelName());
        this.substModel = substModel;
        this.eigenDecomposition = substModel.getEigenDecomposition();
        stateCount = substModel.getDataType().getStateCount();
        markovJumpsCore = new MarkovJumpsCore(stateCount);
        this.type = type;
        setupStorage();
        addModel(substModel);
        dataType = substModel.getDataType();
    }

    protected void setupStorage() {
        rateMatrix = new double[stateCount * stateCount];
        transitionProbs = new double[stateCount * stateCount];
        rateReg = new double[stateCount * stateCount];
        if (PRECOMPUTE) {
            ievcRateRegEvec = new double[stateCount * stateCount];
            tmp1 = new double[stateCount * stateCount];
        }
        registration = new double[stateCount * stateCount];
        reward = new double[stateCount];
    }

    public MarkovJumpsType getType() {
        return type;
    }

    public void setRegistration(double[] inRegistration) {

        if (type == MarkovJumpsType.COUNTS || type == MarkovJumpsType.HISTORY) {

            System.arraycopy(inRegistration, 0, registration, 0, stateCount * stateCount);
            for (int i = 0; i < stateCount; i++) {
                registration[i * stateCount + i] = 0;  // diagonals are zero
            }

        } else if (type == MarkovJumpsType.REWARDS) {

            int index = 0;
            for (int i = 0; i < stateCount; i++) {
                reward[i] = inRegistration[i];
                for (int j = 0; j < stateCount; j++) {
                    if (i == j) {
                        registration[index] = inRegistration[i];
                    } else {
                        registration[index] = 0; // Off-diagonals are zero
                    }
                    index++;
                }
            }

        } else {
            throw new RuntimeException("Unknown expectation type in MarkovJumps");
        }
        regRateChanged = true;
    }

    public double[] getRegistration() {
        return registration;
    }

    private void makeRateRegistrationMatrix(double[] registration,
                                            double[] rateReg,
                                            double[] ievcRateRegEvec) {

        if (type == MarkovJumpsType.COUNTS || type == MarkovJumpsType.HISTORY) {

            substModel.getInfinitesimalMatrix(rateMatrix);
            int index = 0;
            for (int i = 0; i < stateCount; i++) {
                for (int j = 0; j < stateCount; j++) {
                    rateReg[index] = rateMatrix[index] * registration[index];
                    index++;
                }
            }

        } else if (type == MarkovJumpsType.REWARDS) {

            System.arraycopy(registration, 0, rateReg, 0, stateCount * stateCount);

        } else {
            throw new RuntimeException("Unknown expectation type in MarkovJumps");
        }

        if (PRECOMPUTE) {
//            matrixMultiply(rateReg, evec, stateCount, tmp1);
//            matrixMultiply(ievc, tmp1, stateCount, tmp2);
            MarkovJumpsCore.matrixMultiply(rateReg, eigenDecomposition.getEigenVectors(),
                    stateCount, tmp1);
            MarkovJumpsCore.matrixMultiply(eigenDecomposition.getInverseEigenVectors(), tmp1,
                    stateCount, ievcRateRegEvec);
        }

        regRateChanged = false;
    }

    public double getMarginalRate() {

        if (regRateChanged) {
            makeRateRegistrationMatrix(registration, rateReg, ievcRateRegEvec);
        }

        FrequencyModel freqModel = substModel.getFrequencyModel();
        double rate = 0;
        int index = 0;
        for (int i = 0; i < stateCount; i++) {
            double freq_i = freqModel.getFrequency(i);
            for (int j = 0; j < stateCount; j++) {
                rate += freq_i * rateReg[index++];
            }
        }
        return rate;
    }

    public void computeCondStatMarkovJumps(double time,
                                           double[] countMatrix) {

        substModel.getTransitionProbabilities(time, transitionProbs);
        computeCondStatMarkovJumps(time, transitionProbs, countMatrix);
    }


    public double getProcessForSimulant(StateHistory history) {
        final double total;
        if (type == MarkovJumpsType.COUNTS || type == MarkovJumpsType.HISTORY) {
            total = history.getTotalRegisteredCounts(registration);
        } else {
            total = history.getTotalReward(reward);
        }
        return total;
    }


    public void computeCondStatMarkovJumps(double time,
                                           double[] transitionProbs,
                                           double[] countMatrix) {

        if (regRateChanged) {
            makeRateRegistrationMatrix(registration, rateReg, ievcRateRegEvec);
        }

        double[] evec = eigenDecomposition.getEigenVectors();
        double[] ievc = eigenDecomposition.getInverseEigenVectors();
        double[] eval = eigenDecomposition.getEigenValues();

        if (PRECOMPUTE) {
            markovJumpsCore.computeCondStatMarkovJumpsPrecompute(
                    evec, ievc, eval, ievcRateRegEvec, time, transitionProbs, countMatrix);
        } else {
            markovJumpsCore.computeCondStatMarkovJumps(evec, ievc, eval, rateReg, time, transitionProbs, countMatrix);
        }
    }

    public void computeJointStatMarkovJumps(double time,
                                            double[] countMatrix) {

        if (regRateChanged) {
            makeRateRegistrationMatrix(registration, rateReg, ievcRateRegEvec);
        }

        double[] evec = eigenDecomposition.getEigenVectors();
        double[] ievc = eigenDecomposition.getInverseEigenVectors();
        double[] eval = eigenDecomposition.getEigenValues();

        if (PRECOMPUTE) {
            markovJumpsCore.computeJointStatMarkovJumpsPrecompute(evec, ievc, eval, ievcRateRegEvec, time, countMatrix);
        } else {
            markovJumpsCore.computeJointStatMarkovJumps(evec, ievc, eval, rateReg, time, countMatrix);
        }
    }

    public SubstitutionModel getSubstitutionModel() {
        return substModel;
    }

    protected void handleModelChangedEvent(Model model, Object object, int index) {
        if (model == substModel) {
            regRateChanged = true;
        }
    }

    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        // Do nothing
    }

    protected void storeState() {
        // Do nothing
    }

    protected void restoreState() {
        // Do nothing
    }

    protected void acceptState() {
        // Do nothing
    }

    public int stateCount;
    private double[] rateReg;
    private double[] ievcRateRegEvec;
    private double[] tmp1;
    private double[] transitionProbs;
    private double[] rateMatrix;
    protected double[] reward;
    protected double[] registration;

    protected SubstitutionModel substModel;
    private EigenDecomposition eigenDecomposition;
    private MarkovJumpsCore markovJumpsCore;

    private boolean regRateChanged = true;

    protected MarkovJumpsType type;
    protected DataType dataType;

    private static final boolean PRECOMPUTE = true;
}

