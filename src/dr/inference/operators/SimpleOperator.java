/*
 * SimpleOperator.java
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

/**
 *
 */
package dr.inference.operators;

/**
 * @author Sebastian Hoehna
 *
 */
public abstract class SimpleOperator implements MCMCOperator {

    /**
     *
     */
    public SimpleOperator() {
        // TODO Auto-generated constructor stub
    }

    public abstract String getOperatorName();

    /**
     * @return the weight of this operator.
     */
    public final double getWeight() {
        return weight;
    }

    /**
     * Sets the weight of this operator.
     */
    public final void setWeight(double w) {
        if (w > 0) {
            weight = w;
        } else throw new IllegalArgumentException("Weight must be a positive real. (called with " + w + ")");
    }

    public void accept(double deviation) {
        lastDeviation = deviation;

        if (!operateAllowed) {
            operateAllowed = true;
            accepted += 1;
            sumDeviation += deviation;
        } else throw new RuntimeException("Accept/reject methods called twice without operate called in between!");
    }

    public void reject() {
        if (!operateAllowed) {
            operateAllowed = true;
            rejected += 1;
        } else throw new RuntimeException("Accept/reject methods called twice without operate called in between!");
    }

    public final void reset() {
        operateAllowed = true;
        accepted = 0;
        rejected = 0;
        lastDeviation = 0.0;
        sumDeviation = 0.0;
    }

    public final long getCount() {
        return accepted + rejected;
    }

    public final long getAcceptCount() {
        return accepted;
    }

    public final void setAcceptCount(long accepted) {
        this.accepted = accepted;
    }

    public final long getRejectCount() {
        return rejected;
    }

    public final void setRejectCount(long rejected) {
        this.rejected = rejected;
    }

    public final double getMeanDeviation() {
        return sumDeviation / accepted;
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

    public final double getAcceptanceProbability() {
        return (double) accepted / (double) (accepted + rejected);
    }

    public final double getSmoothedAcceptanceProbability() {
        throw new UnsupportedOperationException("Not implemented for SimpleOperator");
    }

    public double getSpan(boolean reset) {
        double span = 0;
        if( spanDeviation[1] > spanDeviation[0] && spanCount > 20 ) {
            span = spanDeviation[1] - spanDeviation[0];

            if( reset ) {
                spanDeviation[0] = Double.MAX_VALUE;
                spanDeviation[1] = -Double.MAX_VALUE;
                spanCount = 0;
            }
        }
        return span;
    }

    public double getMeanEvaluationTime() {
        return (double) sumEvaluationTime / (double) (accepted + rejected);
    }

    public long getTotalEvaluationTime() {
        return sumEvaluationTime;
    }

    public void addEvaluationTime(long time) {
        sumEvaluationTime += time;
    }

    private long sumEvaluationTime = 0;

    private double[] spanDeviation = {Double.MAX_VALUE,-Double.MAX_VALUE};
    private int spanCount = 0;

    private double weight = 1.0;
    private long accepted = 0;
    private long rejected = 0;
    private double sumDeviation = 0.0;
    private double lastDeviation = 0.0;
    protected boolean operateAllowed = true;

}
