/*
 * TwoPhaseOperator.java
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

package dr.inference.operators;

import java.util.ArrayList;
import java.util.List;

import dr.inference.model.Parameter;

/**
 * This class allows to use two different sets of operators.
 * 
 * @author Guy Baele
 */

public class TwoPhaseOperator extends AbstractCoercableOperator {

    public static final boolean DEBUG = false;
    public static final boolean PROVIDE_SAMPLES = false;

    /*private AbstractCoercableOperator[] phaseOneOperators;
    private AbstractCoercableOperator[] phaseTwoOperators;

    double phaseOneTotalWeight, phaseTwoTotalWeight;*/

    private List<AbstractCoercableOperator> phaseOneOperators;
    private List<AdaptableVarianceMultivariateNormalOperator> phaseTwoOperators;
    //private List<AbstractCoercableOperator> phaseTwoOperators;
    private List<AbstractCoercableOperator> currentOperators;

    private SimpleOperatorSchedule phaseOneScheduler;
    private SimpleOperatorSchedule phaseTwoScheduler;
    private SimpleOperatorSchedule currentOperatorScheduler;

    private List<Parameter> parameters;
    private ArrayList<ArrayList<Double>> storedValues;

    private int initial;
    private int burnin;
    private int numberOfCalls;
    private int currentOperatorIndex;

    private boolean switchOperators;

    public TwoPhaseOperator(List<AbstractCoercableOperator> phaseOneOperators, List<AdaptableVarianceMultivariateNormalOperator> phaseTwoOperators, List<Parameter> parameters, int initial, int burnin, double weight, CoercionMode mode) {

        super(mode);

        /*this.phaseOneOperators = phaseOneOperators;
        this.phaseTwoOperators = phaseTwoOperators;

        this.phaseOneTotalWeight = 0.0;
        for (int i = 0; i < phaseOneOperators.length; i++) {
            this.phaseOneTotalWeight += phaseOneOperators[i].getWeight();
        }
        this.phaseTwoTotalWeight = 0.0;
        for (int i = 0; i < phaseTwoOperators.length; i++) {
            this.phaseTwoTotalWeight += phaseTwoOperators[i].getWeight();
        }*/

        if (DEBUG) {
            System.err.println("\nConstructing TwoPhaseOperator");
        }

        setWeight(weight);

        this.initial = initial;
        this.burnin = burnin;
        this.numberOfCalls = 0;

        this.switchOperators = false;

        this.phaseOneOperators = phaseOneOperators;
        this.phaseTwoOperators = phaseTwoOperators;

        this.parameters = parameters;
        this.storedValues = new ArrayList<ArrayList<Double>>();
        for (int i = 0; i < phaseOneOperators.size(); i++) {
            this.storedValues.add(new ArrayList<Double>());
        }

        phaseOneScheduler = new SimpleOperatorSchedule();
        for (MCMCOperator operator : phaseOneOperators) {
            phaseOneScheduler.addOperator(operator);
        }

        if (DEBUG) {
            System.err.println("Phase One Scheduler initiated with size: " + phaseOneScheduler.getOperatorCount()); 
        }

        phaseTwoScheduler = new SimpleOperatorSchedule();
        for (MCMCOperator operator : phaseTwoOperators) {
            phaseTwoScheduler.addOperator(operator);
        }

        if (DEBUG) {
            System.err.println("Phase Two Scheduler initiated with size: " + phaseTwoScheduler.getOperatorCount()); 
        }

        currentOperatorScheduler = phaseOneScheduler;
        currentOperators = phaseOneOperators;

    }

    public double doOperation() {

        if (DEBUG) {
            System.err.println("\nTwoPhaseOperator: doOperation() called");
        }

        numberOfCalls++;
        if (DEBUG) {
            System.err.println("Number of times called: " + numberOfCalls);
        }

        currentOperatorIndex = currentOperatorScheduler.getNextOperatorIndex();
        if (DEBUG) {
            System.err.println("current operator index: " + currentOperatorIndex);
        }

        //don't store anything in the first set of operators themselves
        //store everything in this class to not clutter AVMVN operator with excessive code
        if (numberOfCalls > burnin && !switchOperators) {
            //assume a 1-on-1 relationship between the parameter list and the first phase of operators
            //i.e. each parameter has 1 operator acting on it
            //now we can use currentOperatorIndex to help with the bookkeeping

            //first decide to which of the phase two operators the parameter value needs to be written to
            int phaseTwoCounter = 0;

            //at the same time decide where it actually came from in order to determine its actual value
            //i.e. Parameter might be a CompoundParameter, which complicates things
            int parameterIndex = currentOperatorIndex;
            for (int i = 0; i < phaseTwoOperators.size(); i++) {
                //TODO: this may rely on the AVMVN operator only having 1 CompoundParameter
                if (currentOperatorIndex < phaseTwoOperators.get(i).getParameter().getSize()) {
                    break;
                } else {
                    parameterIndex -= phaseTwoOperators.get(i).getParameter().getSize();
                    phaseTwoCounter++;
                }
            }
            storedValues.get(currentOperatorIndex).add(parameters.get(phaseTwoCounter).getParameterValue(parameterIndex));

            if (DEBUG) {
                System.err.println("Storing values in TwoPhaseOperator");
                System.err.println("currentOperatorIndex: " + currentOperatorIndex);
                System.err.println("parameterIndex: " + parameterIndex);
                System.err.print("storage dimensions: " + storedValues.size());
                for (int i = 0; i < storedValues.size(); i++) {
                    System.err.print(" -> " + storedValues.get(i).size());
                }
                System.err.println();
            }

            /*if (DEBUG) {
                System.err.println("Passing values to phase two operator(s)");
                System.err.println("currentOperatorIndex: " + currentOperatorIndex);
                System.err.println("AVMVN operator assigned: " + phaseTwoCounter);
                System.err.println("parameterIndex: " + parameterIndex);
            }
            phaseTwoOperators.get(phaseTwoCounter).setSample(parameterIndex, parameters.get(phaseTwoCounter).getParameterValue(parameterIndex));
             */
        }

        if (numberOfCalls > initial && !switchOperators) {
            if (DEBUG) {
                System.err.println("Switch from phase one scheduler to phase two scheduler");
            }
            currentOperatorScheduler = phaseTwoScheduler;
            //TODO: fix Java type safety problem below
            currentOperators = (List<AbstractCoercableOperator>)(List<?>) phaseTwoOperators;
            //an extra draw is needed here
            currentOperatorIndex = currentOperatorScheduler.getNextOperatorIndex();

            if (PROVIDE_SAMPLES) {
                //call methods to calculate means and covariance matrix and pass them on to AVMVN operator(s)
                //need to create the appropriate list of lists to pass on to AVMVN operator
                for (int i = 0; i < phaseTwoOperators.size(); i++) {

                    int listSize = phaseTwoOperators.get(i).getParameter().getDimension();

                    ArrayList<ArrayList<Double>> temp = new ArrayList<ArrayList<Double>>();
                    for (int j = 0; j < listSize; j++) {
                        temp.add(new ArrayList<Double>());
                        temp.set(j, storedValues.get(i*phaseTwoOperators.get(i).getParameter().getSize()+j));
                    }

                    phaseTwoOperators.get(i).provideSamples(temp);
                }
            }

            switchOperators = true;
        }

        double logJacobian = (currentOperators.get(currentOperatorIndex)).doOperation();

        return logJacobian;
    }

    /*@Override
    public void accept(double deviation) {
        if (DEBUG) {
            System.err.println("TwoPhaseOperator: accept(double deviation) called");
        }
        currentOperators.get(currentOperatorIndex).accept(deviation);
    }

    @Override
    public void reject() {
        if (DEBUG) {
            System.err.println("TwoPhaseOperator: reject() called");
        }
        currentOperators.get(currentOperatorIndex).reject();
    }

    @Override
    public void reset() {
        if (DEBUG) {
            System.err.println("TwoPhaseOperator: reset() called");
        }
        currentOperators.get(currentOperatorIndex).reset();
    }*/

    public double getCoercableParameter() {
        return currentOperators.get(currentOperatorIndex).getCoercableParameter();
    }

    public void setCoercableParameter(double value) {
        currentOperators.get(currentOperatorIndex).setCoercableParameter(value);
    }

    public double getRawParameter() {
        return currentOperators.get(currentOperatorIndex).getRawParameter();
    }

    public String getPerformanceSuggestion() {
        return currentOperators.get(currentOperatorIndex).getPerformanceSuggestion();
    }

    @Override
    public double getTargetAcceptanceProbability() {
        return currentOperators.get(currentOperatorIndex).getTargetAcceptanceProbability();
    }

    public String getOperatorName() {
        return "twoPhaseOperator(use at own risk)";
    }

}
