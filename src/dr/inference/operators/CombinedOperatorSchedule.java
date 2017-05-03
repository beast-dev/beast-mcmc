/*
 * CombinedOperatorSchedule.java
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

import dr.math.MathUtils;

import java.util.List;
import java.util.Vector;

/**
 * Package: CombinedOperatorSchedule
 * Description:
 * <p/>
 * <p/>
 * Created by
 *
 * @author Alexander V. Alekseyenko (alexander.alekseyenko@gmail.com)
 *         Date: Dec 23, 2009
 *         Time: 10:17:28 AM
 */
public class CombinedOperatorSchedule implements OperatorSchedule {

    List<OperatorSchedule> operatorSchedules = null;

    public CombinedOperatorSchedule() {
        operatorSchedules = new Vector<OperatorSchedule>();
    }

    public void addOperatorSchedule(OperatorSchedule os) {
        operatorSchedules.add(os);
    }

    public int getScheduleCount() {
        return operatorSchedules.size();
    }

    public int getNextOperatorIndex() {

        final int v = MathUtils.nextInt(operatorSchedules.size());

        int offset = 0;
        for (int i = 0; i < v; ++i) {
            offset += operatorSchedules.get(i).getOperatorCount();
        }

        return offset + operatorSchedules.get(v).getNextOperatorIndex();
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
        System.err.println("ERROR: addOperator() should not be called on CombinedOperatorSchedule, " +
                "use individual operator schedules to add operators. This call will have no affect!");
    }

    public void addOperators(List<MCMCOperator> v) {
        System.err.println("ERROR: addOperators() should not be called on CombinedOperatorSchedule, " +
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
        long minCount = Long.MAX_VALUE;
        for (OperatorSchedule os : operatorSchedules) {
            if (os.getMinimumAcceptAndRejectCount() < minCount) {
                minCount = os.getMinimumAcceptAndRejectCount();
            }
        }
        return minCount;
    }


}
