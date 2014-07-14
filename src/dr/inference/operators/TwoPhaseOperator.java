/*
 * TwoPhaseOperator.java
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

package dr.inference.operators;

import java.util.List;

/**
 * This class allows to use two different sets of operators.
 * 
 * @author Guy Baele
 */

public class TwoPhaseOperator extends AbstractCoercableOperator {

    public static final boolean DEBUG = true;
    
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
    
    private int initial;
    private int burnin;
    private int numberOfCalls;
    private int currentOperatorIndex;
    
    public TwoPhaseOperator(List<AbstractCoercableOperator> phaseOneOperators, List<AdaptableVarianceMultivariateNormalOperator> phaseTwoOperators, int initial, int burnin, double weight, CoercionMode mode) {
        
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
        
        this.phaseOneOperators = phaseOneOperators;
        this.phaseTwoOperators = phaseTwoOperators;
        
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
    
    public double doOperation() throws OperatorFailedException {
        
        if (DEBUG) {
            System.err.println("\nTwoPhaseOperator: doOperation() called");
        }
        
        numberOfCalls++;
        if (DEBUG) {
            System.err.println("Number of times called: " + numberOfCalls);
        }
        
        //SHOULD THIS BE PLACED HERE? WHEN IN DOOPERATION() DOES AVMVN STORE VALUES?
        /*if (numberOfCalls > burnin) {
            if (DEBUG) {
                System.err.println("Start passing values to phase two operator(s) if operator(s) allow(s) this");
            }
        }*/
        
        if (numberOfCalls > initial) {
            if (DEBUG) {
                System.err.println("Initialize AVMVN operator(s) with stored samples");
            }
            //copy stored samples to the AVMVN operator(s)
            
            
            
            
            
            if (DEBUG) {
                System.err.println("Switch from phase one scheduler to phase two scheduler");
            }
            currentOperatorScheduler = phaseTwoScheduler;
            //TODO: fix Java type safety problem below
            currentOperators = (List<AbstractCoercableOperator>)(List<?>) phaseTwoOperators;
        }
        
        currentOperatorIndex = currentOperatorScheduler.getNextOperatorIndex();
        
        if (DEBUG) {
            System.err.println("current operator index: " + currentOperatorIndex);
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
