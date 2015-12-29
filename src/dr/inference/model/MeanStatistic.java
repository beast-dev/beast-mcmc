/*
 * MeanStatistic.java
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

import java.util.ArrayList;
import java.util.List;

/**
 * @version $Id: MeanStatistic.java,v 1.9 2005/05/24 20:26:00 rambaut Exp $
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 */
public class MeanStatistic extends Statistic.Abstract {

	public MeanStatistic(String name) {
		super(name);
	}

	public void addStatistic(Statistic statistic) {
		statistics.add(statistic);
	}

	public int getDimension() { return 1; }

	/** @return mean of contained statistics */
	public double getStatisticValue(int dim) {
		double sum = 0.0;
		int dimensionCount = 0;
		int n;

		for (Statistic statistic : statistics) {
			n = statistic.getDimension();

			for (int j = 0; j < n; j++) {
				sum += statistic.getStatisticValue(j);
			}
			dimensionCount += n;
		}

		return sum / dimensionCount;
	}

	// ****************************************************************
	// Private and protected stuff
	// ****************************************************************

	private List<Statistic> statistics = new ArrayList<Statistic>();
}
