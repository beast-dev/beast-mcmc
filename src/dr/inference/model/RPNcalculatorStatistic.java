/*
 * RPNcalculatorStatistic.java
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

package dr.inference.model;

import java.util.Map;

/**
 *  A statistic based on evaluating simple expressions.
 *
 * The expressions are in RPN, so no parsing issues. whitspace separated. Variables (other statistics),
 * constants and operations. Currently just the basic four, but easy to extend.
 *
 * @author Joseph Heled
 */
public class RPNcalculatorStatistic extends Statistic.Abstract {

    private final RPNexpressionCalculator[] expressions;
    private final String[] names;
    private final Map<String, Statistic> variables;

    RPNexpressionCalculator.GetVariable vars = new RPNexpressionCalculator.GetVariable() {
        public double get(String name) {
            return variables.get(name).getStatisticValue(0);
        }
    };

    public RPNcalculatorStatistic(String name, String[] expressions, String[] names,
                                  Map<String, Statistic> variables) {
        super(name);

        this.expressions = new RPNexpressionCalculator[expressions.length];
        for(int i = 0; i < expressions.length; ++ i) {
            this.expressions[i] = new RPNexpressionCalculator(expressions[i]);

            String err = this.expressions[i].validate();
            if( err != null ) {
                throw new RuntimeException("Error in expression " + i + ": " + err);
            }
        }

        this.names = names;
        this.variables = variables;
    }

	public int getDimension() {
        return expressions.length;
    }

    public String getDimensionName(int dim) {
        return names[dim];
    }

    /** @return the value of the expression */
	public double getStatisticValue(int dim) {
        return expressions[dim].evaluate(vars);
	}

}