/*
 * PrefetchOperatorSchedule.java
 *
 * Copyright (c) 2002-2017 Alexei Drummond, Andrew Rambaut and Marc Suchard
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
 * PrefetchOperatorSchedule
 * @author rambaut
 */
public class PrefetchOperatorSchedule implements OperatorSchedule {

    private final OperatorSchedule operatorSchedule;

    private MCMCOperator currentPrefetchOperator = null;

    public PrefetchOperatorSchedule(OperatorSchedule operatorSchedule) {
        this.operatorSchedule = operatorSchedule;
    }

    public int getNextOperatorIndex() {
        return operatorSchedule.getNextOperatorIndex();
    }

    public int getOperatorCount() {
        return operatorSchedule.getOperatorCount();
    }

    public MCMCOperator getOperator(int index) {
        MCMCOperator operator;

        if (currentPrefetchOperator != null) {
            operator = currentPrefetchOperator;
            if (((Prefetchable)currentPrefetchOperator).prefetchingDone()) {
                currentPrefetchOperator = null;
            }
        } else {
            operator = operatorSchedule.getOperator(index);

            if (operator instanceof Prefetchable) {
                currentPrefetchOperator = operator;
            }
        }

        return operator;
    }

    public void addOperator(MCMCOperator op) {
        System.err.println("ERROR: addOperator() should not be called on PrefetchOperatorSchedule, " +
                "use individual operator schedules to add operators. This call will have no affect!");
    }

    public void addOperators(List<MCMCOperator> v) {
        System.err.println("ERROR: addOperators() should not be called on PrefetchOperatorSchedule, " +
                "use individual operator schedules to add operators. This call will have no affect!");
    }

    public void operatorsHaveBeenUpdated() {
        operatorSchedule.operatorsHaveBeenUpdated();
    }

    public double getOptimizationTransform(double d) {
        return operatorSchedule.getOptimizationTransform(d);
    }

    public int getMinimumAcceptAndRejectCount() {
        return operatorSchedule.getMinimumAcceptAndRejectCount();
    }


}
