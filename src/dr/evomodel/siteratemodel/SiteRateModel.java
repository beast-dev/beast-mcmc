/*
 * SiteRateModel.java
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

package dr.evomodel.siteratemodel;

import dr.inference.model.Model;

/**
 * @author Andrew Rambaut
 * @author Marc A. Suchard
 */

public interface SiteRateModel extends Model {

    /**
     * @return the number of categories of substitution processes
     */
    int getCategoryCount();

    /**
     * Get an array of the relative rates of sites in each category. These
     * may include the 'mu' parameter, an overall scaling of the siteRateModel.
     *
     * @return an array of the rates.
     */
    double[] getCategoryRates();

    /**
     * Get an array of the expected proportion of sites in each category.
     *
     * @return an array of the proportions.
     */
    double[] getCategoryProportions();

    /**
     * Get the rate for a particular category. This may include the 'mu'
     * parameter, an overall scaling of the siteRateModel.
     *
     * @param category the category number
     * @return the rate.
     */
    double getRateForCategory(int category);

    /**
     * Get the expected proportion of sites in this category.
     *
     * @param category the category number
     * @return the proportion.
     */
    double getProportionForCategory(int category);

}
