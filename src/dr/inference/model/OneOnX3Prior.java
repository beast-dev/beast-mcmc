/*
 * OneOnX3Prior.java
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

import java.util.ArrayList;

/**
 * @author Chieh-Hsi Wu
 *
 * This an improper prior 1/(x_i^3) which is the Jeffrey's prior for parameters of Normal distribution
 * (and Log-normal distribution?)
 *
 */
public class OneOnX3Prior extends Likelihood.Abstract{

    public OneOnX3Prior() {

        super(null);
    }

    /**
     * Adds a statistic, this is the data for which the Prod_i (1/x_i^3) prior is calculated.
     *
     * @param data the statistic to compute density of
     */
    public void addData(Statistic data) {
        dataList.add(data);
    }


    protected ArrayList<Statistic> dataList = new ArrayList<Statistic>();

    /**
     * Overridden to always return false.
     */
    protected boolean getLikelihoodKnown() {
        return false;
    }

    /**
     * Calculate the log likelihood of the current state.
     *
     * @return the log likelihood.
     */
    public double calculateLogLikelihood() {

        double logL = 0.0;

        for (Statistic statistic : dataList) {
            for (int j = 0; j < statistic.getDimension(); j++) {
                logL -= 3*Math.log(statistic.getStatisticValue(j));
            }
        }
        return logL;
    }


    public String prettyName() {
        String s = "OneOnX3" + "(";
        for (Statistic statistic : dataList) {
            s = s + statistic.getStatisticName() + ",";
        }
        return s.substring(0, s.length() - 1) + ")";
    }
}
