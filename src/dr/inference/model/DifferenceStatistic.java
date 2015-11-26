/*
 * DifferenceStatistic.java
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
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: ProductStatistic.java,v 1.2 2005/05/24 20:26:00 rambaut Exp $
 */
public class DifferenceStatistic extends Statistic.Abstract {

    private final boolean absolute;

    public DifferenceStatistic(String name, Statistic term1, Statistic term2, boolean absolute) {
        super(name);

        this.term1 = term1;
        this.term2 = term2;

        if (term1.getDimension() != 1 &&
                term2.getDimension() != 1 &&
                term1.getDimension() != term2.getDimension()) {
            throw new IllegalArgumentException();
        }

        if (term2.getDimension() == 1) {
            dimension = term1.getDimension();
        } else {
            dimension = term2.getDimension();
        }

        this.absolute = absolute;
    }

    public int getDimension() {
        return dimension;
    }

    /**
     * @return mean of contained statistics
     */
    public double getStatisticValue(int dim) {

        double statistic;

        if (term1.getDimension() == 1) {
            statistic = term1.getStatisticValue(0) - term2.getStatisticValue(dim);
        } else if (term2.getDimension() == 1) {
            statistic = term1.getStatisticValue(dim) - term1.getStatisticValue(0);
        } else {
            statistic = term1.getStatisticValue(dim) - term2.getStatisticValue(dim);
        }
        if (absolute) statistic = Math.abs(statistic);

        return statistic;
    }

    // ****************************************************************
    // Private and protected stuff
    // ****************************************************************

    private int dimension = 0;
    private Statistic term1 = null;
    private Statistic term2 = null;
}
