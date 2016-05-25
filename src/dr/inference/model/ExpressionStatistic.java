/*
 * ExpressionStatistic.java
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

import java.util.Vector;

/**
 * @version $Id: ExpressionStatistic.java,v 1.4 2005/05/24 20:26:00 rambaut Exp $
 *
 * @author Alexei Drummond
 */
public class ExpressionStatistic extends Statistic.Abstract {

    String expression = "";

	public ExpressionStatistic(String name, String expression) {
		super(name);
        this.expression = expression;
	}
	
	public void addStatistic(Statistic statistic) {
        if (statistic.getDimension() != 1) {
            throw new IllegalArgumentException("Can only have statistics of dimension 1");
        }

		statistics.add(statistic);
	}
	
	public int getDimension() { return 1; }

	/** @return the value of the expression */
	public double getStatisticValue(int dim) {

		System.err.println("Error in parsing expression " + expression + " : JEP expression parser not included with this version");
        return 0;
	}

	// ****************************************************************
	// Private and protected stuff
	// ****************************************************************
	
	private Vector statistics = new Vector();
}
