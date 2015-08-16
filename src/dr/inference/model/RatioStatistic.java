/*
 * RatioStatistic.java
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

/**
 * @version $Id: ProductStatistic.java,v 1.2 2005/05/24 20:26:00 rambaut Exp $
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 */
public class RatioStatistic extends Statistic.Abstract {

	public RatioStatistic(String name, Statistic numerator, Statistic denominator) {
		super(name);

        this.numerator = numerator;
        this.denominator = denominator;

        if (denominator.getDimension() != 1 &&
                numerator.getDimension() != 1 &&
                denominator.getDimension() != numerator.getDimension()) {
            throw new IllegalArgumentException();
        }

        if (denominator.getDimension() == 1) {
            dimension = numerator.getDimension();
        } else {
            dimension = denominator.getDimension();
        }
    }

    public int getDimension() { return dimension; }

	/** @return mean of contained statistics */
	public double getStatisticValue(int dim) {

        if (numerator.getDimension() == 1) {
            return numerator.getStatisticValue(0) / denominator.getStatisticValue(dim);
        } else if (denominator.getDimension() == 1) {
            return numerator.getStatisticValue(dim) / denominator.getStatisticValue(0);
        } else {
            return numerator.getStatisticValue(dim) / denominator.getStatisticValue(dim);
        }
    }

	// ****************************************************************
	// Private and protected stuff
	// ****************************************************************

    private int dimension = 0;
    private Statistic numerator = null;
	private Statistic denominator = null;
}
