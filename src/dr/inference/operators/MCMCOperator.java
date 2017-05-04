/*
 * MCMCOperator.java
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

import java.io.Serializable;

/**
 * An MCMC operator.
 *
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @version $Id: MCMCOperator.java,v 1.6 2005/06/14 10:40:34 rambaut Exp $
 */
public interface MCMCOperator extends Serializable {

    public static final String WEIGHT = "weight";

// This attribute is now called AUTO_OPTIMIZE and is in CoercableMCMCOperator
//	public static final String ADAPT = "adapt";

    /**
     * operates on the model.
     *
     * @return the hastings ratio of this operator.
     */
    double operate();

    /**
     * Called to tell operator that operation was accepted
     *
     * @param deviation the log ratio accepted on
     */
    void accept(double deviation);

    /**
     * Called to tell operator that operation was rejected
     */
    void reject();

    /**
     * Reset operator acceptance records.
     */
    void reset();

    /**
     * @return the total number of operations since last call to reset().
     */
    long getCount();

    /**
     * @return the number of acceptances since last call to reset().
     */
    long getAcceptCount();

    /**
     * Set the number of acceptances since last call to reset(). This is used
     * to restore the state of the operator
     *
     * @param acceptCount number of acceptances
     */
    void setAcceptCount(long acceptCount);

    /**
     * @return the number of rejections since last call to reset().
     */
    long getRejectCount();

    /**
     * Set the number of rejections since last call to reset(). This is used
     * to restore the state of the operator
     *
     * @param rejectCount number of rejections
     */
    void setRejectCount(long rejectCount);

    /**
     * @return the mean deviation in log posterior per accepted operations.
     */
    double getMeanDeviation();

    double getSumDeviation();

    //double getSpan(boolean reset);

    void setSumDeviation(double sumDeviation);

    /**
     * @return the optimal acceptance probability
     */
    double getTargetAcceptanceProbability();

    /**
     * @return the minimum acceptable acceptance probability
     */
    double getMinimumAcceptanceLevel();

    /**
     * @return the maximum acceptable acceptance probability
     */
    double getMaximumAcceptanceLevel();

    /**
     * @return the minimum good acceptance probability
     */
    double getMinimumGoodAcceptanceLevel();

    /**
     * @return the maximum good acceptance probability
     */
    double getMaximumGoodAcceptanceLevel();

    /**
     * @return a short descriptive message of the performance of this operator.
     */
    String getPerformanceSuggestion();

    /**
     * @return the relative weight of this operator.
     */
    double getWeight();

    /**
     * sets the weight of this operator. The weight
     * determines the proportion of time spent using
     * this operator. This is relative to a 'standard'
     * operator weight of 1.
     *
     * @param weight the relative weight of this parameter - should be positive.
     */
    void setWeight(double weight);

    /**
     * @return the name of this operator
     */
    String getOperatorName();

    double getMeanEvaluationTime();

    void addEvaluationTime(long time);

    long getTotalEvaluationTime();

    class Utils {

        public static double getAcceptanceProbability(MCMCOperator op) {
            final long accepted = op.getAcceptCount();
            final long rejected = op.getRejectCount();
            return (double) accepted / (double) (accepted + rejected);
        }

        public static long getOperationCount(MCMCOperator op) {
            return op.getAcceptCount() + op.getRejectCount();
        }
    }
}
