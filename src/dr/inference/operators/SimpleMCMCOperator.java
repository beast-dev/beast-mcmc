/*
 * SimpleMCMCOperator.java
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

import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;

import java.util.ArrayDeque;
import java.util.Deque;

public abstract class SimpleMCMCOperator implements MCMCOperator {

    private final static int SMOOTHED_ACCEPTANCE_WINDOW_SIZE = 100;

    public abstract String getOperatorName();

    public SimpleMCMCOperator() {
    }

    public SimpleMCMCOperator(double weight) {
        this.weight = weight;
    }

    /**
     * @return the weight of this operator.
     */
    public final double getWeight() {
        return weight;
    }

    public void setPathParameter(double beta) {
        throw new IllegalArgumentException("Path parameter has no effect on Metropolis-Hastings kernels." +
                "\nGibbs samplers need an implementation for use in power-posteriors");
    }

    /**
     * Sets the weight of this operator.
     */
    public final void setWeight(double w) {
        if( w > 0 ) {
            weight = w;
        } else {
            throw new IllegalArgumentException(
                    "Weight must be a positive real, but tried to set weight to "
                            + w);
        }
    }

    public void accept(double deviation) {
        lastDeviation = deviation;

        if( !operateAllowed ) {
            operateAllowed = true;
            acceptCount += 1;
            sumDeviation += deviation;

            windowAcceptance.addLast(1);
            if (windowAcceptance.size() > SMOOTHED_ACCEPTANCE_WINDOW_SIZE) {
                windowAcceptance.removeFirst();
            }
        } else {
            throw new RuntimeException(
                    "Accept/reject methods called twice without operate called in between!");
        }
    }


    public void reject() {
        if( !operateAllowed ) {
            operateAllowed = true;
            rejectCount += 1;

            windowAcceptance.addLast(0);
            if (windowAcceptance.size() > SMOOTHED_ACCEPTANCE_WINDOW_SIZE) {
                windowAcceptance.removeFirst();
            }
        } else {
            throw new RuntimeException(
                    "Accept/reject methods called twice without operate called in between!");
        }
    }

    public void reset() {
        operateAllowed = true;
        acceptCount = 0;
        rejectCount = 0;
        lastDeviation = 0.0;
        sumDeviation = 0.0;

        windowAcceptance.clear();
    }

    public final long getCount() {
        return acceptCount + rejectCount;
    }

    public final long getAcceptCount() {
        return acceptCount;
    }

    public final void setAcceptCount(long acceptCount) {
        this.acceptCount = acceptCount;
    }

    public final long getRejectCount() {
        return rejectCount;
    }

    public final void setRejectCount(long rejectCount) {
        this.rejectCount = rejectCount;
    }

    public final double getMeanDeviation() {
        return sumDeviation / acceptCount;
    }

    public final double getDeviation() {
        return lastDeviation;
    }

    public final double getSumDeviation() {
        return sumDeviation;
    }

    public final void setSumDeviation(double sumDeviation) {
        this.sumDeviation = sumDeviation;
    }

    public final double operate() {
        if( operateAllowed ) {
            operateAllowed = false;
            return doOperation();
        } else {
            throw new RuntimeException(
                    "Operate called twice without accept/reject in between!");
        }
    }

    public final double operate(Likelihood likelihood) {
        if( operateAllowed ) {
            operateAllowed = false;
            return doOperation(likelihood);
        } else {
            throw new RuntimeException(
                    "Operate called twice without accept/reject in between!");
        }
    }

    @Override
    public final double getAcceptanceProbability() {
        return (double) acceptCount / (double) (acceptCount + rejectCount);
    }

    @Override
    public final double getSmoothedAcceptanceProbability() {
        int prob = 0;
        for (int accept : windowAcceptance) {
            prob += accept;
        }
        return (double) prob / (double) (windowAcceptance.size());
    }

    /**
     * Called by operate(), does the actual operation.
     *
     * @return the hastings ratio
     */
    public double doOperation(Likelihood likelihood) {
        return 0.0;
    }

    @Override
    public double getMeanEvaluationTime() {
        return (double) sumEvaluationTime / (double) (acceptCount + rejectCount);
    }

    @Override
    public long getTotalEvaluationTime() {
        return sumEvaluationTime;
    }

    @Override
    public void addEvaluationTime(long time) {
        sumEvaluationTime += time;
    }

    @Override
    public double getMeanCalculationCount() {
        return (double) sumCalculationCount / (double) (acceptCount + rejectCount);
    }

    @Override
    public void addCalculationCount(long count) {
        sumCalculationCount += count;
    }

    @Override
    public long getTotalCalculationCount() {
        return sumCalculationCount;
    }

    /**
     * Called by operate(), does the actual operation.
     *
     * @return the hastings ratio
     */
    public abstract double doOperation();

    public Parameter getParameter() {
        throw new RuntimeException("not implemented for operator of class " + getOperatorName());
    }

    private double weight = 1.0;
    private long acceptCount = 0;
    private long rejectCount = 0;

    private double sumDeviation = 0.0;
    private double lastDeviation = 0.0;

    private boolean operateAllowed = true;

    private long sumEvaluationTime = 0;
    private long sumCalculationCount = 0;

    private Deque<Integer> windowAcceptance = new ArrayDeque<>();
}
