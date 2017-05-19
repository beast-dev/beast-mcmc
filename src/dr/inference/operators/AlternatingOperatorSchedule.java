/*
 * ComplexOperatorSchedule.java
 *
 * Copyright (c) 2002-2016 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

import dr.math.MathUtils;

import java.util.ArrayList;
import java.util.List;
import java.util.Vector;

/**
 * File: AlternatingOperatorSchedule
 * @author Andrew Rambaut
 */
public class AlternatingOperatorSchedule implements OperatorSchedule {

    private final List<OperatorSchedule> operatorSchedules = new ArrayList<OperatorSchedule>();
    private final List<Long> operatorCounts = new ArrayList<Long>();

    private int currentSchedule;
    private long scheduleCount;

    public AlternatingOperatorSchedule() {
        currentSchedule = 0;
        scheduleCount = 0;
    }

    public void addOperatorSchedule(OperatorSchedule os, long count) {
        operatorSchedules.add(os);
        operatorCounts.add(count);
    }

    public int getNextOperatorIndex() {

        if (scheduleCount >= operatorCounts.get(currentSchedule)) {
            currentSchedule += 1;
            if (currentSchedule >= operatorSchedules.size()) {
                currentSchedule = 0;
            }
            scheduleCount = 0;
        }
        scheduleCount += 1;

        int offset = 0;
        for (int i = 0; i < currentSchedule; ++i) {
            offset += operatorSchedules.get(i).getOperatorCount();
        }

        return offset + operatorSchedules.get(currentSchedule).getNextOperatorIndex();
    }

    public void reset() {
        for (OperatorSchedule os : operatorSchedules) {
            for (int i = 0; i < os.getOperatorCount(); ++i) {
                os.getOperator(i).reset();
            }
        }
    }

    public int getOperatorCount() {
        int operatorCount = 0;
        for (OperatorSchedule os : operatorSchedules) {
            operatorCount += os.getOperatorCount();
        }
        return operatorCount;
    }

    public MCMCOperator getOperator(int index) {
        for (OperatorSchedule os : operatorSchedules) {
            int opCount = os.getOperatorCount();
            if (index < opCount) {
                return os.getOperator(index);
            } else {
                index -= opCount;
            }
        }
        // if we reach here the index must be out of bounds return null
        return null;
    }

    public void addOperator(MCMCOperator op) {
        System.err.println("ERROR: addOperator() should not be called on AlternatingOperatorSchedule, " +
                "use individual operator schedules to add operators. This call will have no affect!");
    }

    public void addOperators(List<MCMCOperator> v) {
        System.err.println("ERROR: addOperators() should not be called on AlternatingOperatorSchedule, " +
                "use individual operator schedules to add operators. This call will have no affect!");
    }

    public void operatorsHasBeenUpdated() {
        for (OperatorSchedule os : operatorSchedules) {
            os.operatorsHasBeenUpdated();
        }
    }

    public double getOptimizationTransform(double d) {
        if (operatorSchedules.size() > 0)
            return operatorSchedules.get(0).getOptimizationTransform(d);
        else
            return 0;
    }

    public long getMinimumAcceptAndRejectCount() {
        long minCount = Integer.MAX_VALUE;
        for (OperatorSchedule os : operatorSchedules) {
            if (os.getMinimumAcceptAndRejectCount() < minCount) {
                minCount = os.getMinimumAcceptAndRejectCount();
            }
        }
        return minCount;
    }


}
