/*
 * VarianceStatistic.java
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

import dr.stats.DiscreteStatistics;

import java.util.Vector;

/**
 * @version $Id: VarianceStatistic.java,v 1.3 2005/05/24 20:26:00 rambaut Exp $
 *
 * @author Alexei Drummond
 */
public class VarianceStatistic extends Statistic.Abstract {

	public VarianceStatistic(String name) {
		super(name);
	}
	
	public void addStatistic(Statistic statistic) {
		statistics.add(statistic);
        int dimensionCount = 0;
        for (int i = 0; i < statistics.size(); i++) {
			statistic = (Statistic)statistics.get(i);
			dimensionCount += statistic.getDimension();
        }
        values = new double[dimensionCount];
	}
	
	public int getDimension() { return 1; }

	/** @return variance of contained statistics */
	public final double getStatisticValue(int dim) {
		int n;
		Statistic statistic;

        int index = 0;
        for (int i = 0; i < statistics.size(); i++) {
            statistic = (Statistic)statistics.get(i);
			n = statistic.getDimension();
            for (int j = 0; j < n; j++) {
				values[index] = statistic.getStatisticValue(j);
			    index += 1;
            }
		}
		
		return DiscreteStatistics.variance(values);
	}

	// ****************************************************************
	// Private and protected stuff
	// ****************************************************************
	
	private Vector statistics = new Vector();
    private double[] values = null;
}
