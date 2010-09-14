/*
 * SimpleOperatorSchedule.java
 *
 * Copyright (C) 2002-2006 Alexei Drummond and Andrew Rambaut
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

import dr.inference.loggers.LogColumn;
import dr.inference.loggers.Loggable;
import dr.inference.loggers.NumberColumn;
import dr.math.MathUtils;

import java.util.List;
import java.util.Vector;

/**
 * This class implements a simple operator schedule.
 *
 * @author Alexei Drummond
 * @version $Id: SimpleOperatorSchedule.java,v 1.5 2005/06/14 10:40:34 rambaut Exp $
 */
public class SimpleOperatorSchedule implements OperatorSchedule, Loggable {

	List<MCMCOperator> operators = null;
	double totalWeight = 0;
	int current = 0;
	boolean sequential = false;
	int optimizationSchedule = OperatorSchedule.DEFAULT_SCHEDULE;

	public SimpleOperatorSchedule() {
		operators = new Vector<MCMCOperator>();
	}

	public void addOperators(List<MCMCOperator> operators) {
		for (MCMCOperator operator : operators) {
			operators.add(operator);
			totalWeight += operator.getWeight();
		}
	}

	public void operatorsHasBeenUpdated() {
		totalWeight = 0.0;
		for (MCMCOperator operator : operators) {
			totalWeight += operator.getWeight();
		}
	}

	public void addOperator(MCMCOperator op) {
		operators.add(op);
		totalWeight += op.getWeight();
	}

	public double getWeight(int index) {
		return operators.get(index).getWeight();
	}

	public int getNextOperatorIndex() {

		if (sequential) {
			int index = getWeightedOperatorIndex(current);
			current += 1;
			if (current >= totalWeight) {
				current = 0;
			}
			return index;
		}

        final double v = MathUtils.nextDouble();
        //System.err.println("v=" + v);
        return getWeightedOperatorIndex(v * totalWeight);
	}

	public void setSequential(boolean seq) {
		sequential = seq;
	}

	private int getWeightedOperatorIndex(double q) {
		int index = 0;
		double weight = getWeight(index);
		while (weight <= q) {
			index += 1;
			weight += getWeight(index);
		}
		return index;
	}

	public MCMCOperator getOperator(int index) {
		return operators.get(index);
	}

	public int getOperatorCount() {
		return operators.size();
	}

	public double getOptimizationTransform(double d) {
        switch( optimizationSchedule ) {
            case LOG_SCHEDULE:  return Math.log(d);
            case SQRT_SCHEDULE: return Math.sqrt(d);
        }
		return d;
	}

	public void setOptimizationSchedule(int schedule) {
		optimizationSchedule = schedule;
	}

    public int getMinimumAcceptAndRejectCount() {
        int minCount = Integer.MAX_VALUE;
        for( MCMCOperator op : operators ) {
            if( op.getAcceptCount() < minCount || op.getRejectCount() < minCount ) {
                minCount = op.getCount();
            }
        }
        return minCount;
    }

	// **************************************************************
	// Loggable IMPLEMENTATION
	// **************************************************************

	/**
	 * @return the log columns.
	 */
	public LogColumn[] getColumns() {
		LogColumn[] columns = new LogColumn[getOperatorCount()];
		for (int i = 0; i < getOperatorCount(); i++) {
			MCMCOperator op = getOperator(i);
			columns[i] = new OperatorColumn(op.getOperatorName(), op);
		}
		return columns;
	}

    private class OperatorColumn extends NumberColumn {
		private final MCMCOperator op;

		public OperatorColumn(String label, MCMCOperator op) {
			super(label);
			this.op = op;
		}

		public double getDoubleValue() {
			return MCMCOperator.Utils.getAcceptanceProbability(op);
		}
	}
}
