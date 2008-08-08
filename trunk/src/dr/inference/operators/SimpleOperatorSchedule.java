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
import dr.xml.*;

import java.util.List;
import java.util.Vector;
import java.util.logging.Logger;

/**
 * This class implements a simple operator schedule.
 *
 * @author Alexei Drummond
 * @version $Id: SimpleOperatorSchedule.java,v 1.5 2005/06/14 10:40:34 rambaut Exp $
 */
public class SimpleOperatorSchedule implements OperatorSchedule, Loggable {

	public static final String OPERATOR_SCHEDULE = "operators";
	public static final String SEQUENTIAL = "sequential";
	public static final String OPTIMIZATION_SCHEDULE = "optimizationSchedule";

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

		return getWeightedOperatorIndex(MathUtils.nextDouble() * totalWeight);
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
		if (optimizationSchedule == OperatorSchedule.LOG_SCHEDULE)
			return Math.log(d);
		if (optimizationSchedule == OperatorSchedule.SQRT_SCHEDULE)
			return Math.sqrt(d);
		return d;
	}

	public void setOptimizationSchedule(int schedule) {
		optimizationSchedule = schedule;
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
		private MCMCOperator op;

		public OperatorColumn(String label, MCMCOperator op) {
			super(label);
			this.op = op;
		}

		public double getDoubleValue() {
			return MCMCOperator.Utils.getAcceptanceProbability(op);
		}
	}

	public static dr.xml.XMLObjectParser PARSER = new AbstractXMLObjectParser() {

		public String getParserName() {
			return OPERATOR_SCHEDULE;
		}

		public Object parseXMLObject(XMLObject xo) throws XMLParseException {

			SimpleOperatorSchedule schedule = new SimpleOperatorSchedule();

			if (xo.hasAttribute(SEQUENTIAL)) {
				schedule.setSequential(xo.getBooleanAttribute(SEQUENTIAL));
			}


			if (xo.hasAttribute(OPTIMIZATION_SCHEDULE)) {
				String type = xo.getStringAttribute(OPTIMIZATION_SCHEDULE);
                Logger.getLogger("dr.inference").info("Optimization Scheduule: " + type);
                 
                if (type.equals(OperatorSchedule.LOG_STRING))
					schedule.setOptimizationSchedule(OperatorSchedule.LOG_SCHEDULE);
				else if (type.equals(OperatorSchedule.SQRT_STRING))
					schedule.setOptimizationSchedule(SQRT_SCHEDULE);
				else if (!type.equals(OperatorSchedule.DEFAULT_STRING))
					throw new RuntimeException("Unsupported optimization schedule");
			}

			for (int i = 0; i < xo.getChildCount(); i++) {
				Object child = xo.getChild(i);
				if (child instanceof MCMCOperator) {
					schedule.addOperator((MCMCOperator) child);
				}
			}
			return schedule;
		}

		//************************************************************************
		// AbstractXMLObjectParser implementation
		//************************************************************************

		public XMLSyntaxRule[] getSyntaxRules() {
			return rules;
		}

		private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
				AttributeRule.newBooleanRule(SEQUENTIAL, true),
				new ElementRule(MCMCOperator.class, 1, Integer.MAX_VALUE),
				AttributeRule.newStringRule(OPTIMIZATION_SCHEDULE, true)
		};

		public String getParserDescription() {
			return "A simple operator scheduler";
		}

		public Class getReturnType() {
			return SimpleOperatorSchedule.class;
		}
	};
}
