/*
 * LogarithmStatistic.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
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
 *
 */

package dr.inference.model;

/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 */
public class LogarithmStatistic extends Statistic.Abstract {

    private final Statistic statistic;
    private final double base;

    public LogarithmStatistic(String name, Statistic statistic, double base) {
        super(name);
        this.statistic = statistic;
        this.base = base;
    }

    public int getDimension() {
        return statistic.getDimension();
    }

    /**
     * @return mean of contained statistics
     */
    public double getStatisticValue(int dim) {
        if (base <= 1.0) {
            return Math.log(statistic.getStatisticValue(dim));
        } else {
            return Math.log(statistic.getStatisticValue(dim)) / Math.log(base);
        }
    }

}