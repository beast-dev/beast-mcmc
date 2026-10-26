/*
 * DifferentiableSiteRateModel.java
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

package dr.evomodel.siteratemodel;

import dr.inference.model.Parameter;

/**
 * A site-rate model whose category rates and proportions can be differentiated with respect to its parameters.
 *
 * @author Xiang Ji
 * @author Marc A. Suchard
 */
public interface DifferentiableSiteRateModel extends SiteRateModel {

    /**
     * @return true if getCategoryDifferentials accepts the (scalar) parameter
     */
    boolean isDifferentiableWrtParameter(Parameter parameter);

    /**
     * Fills d rate_c / d parameter and d proportion_c / d parameter, in the category order of getCategoryRates
     */
    void getCategoryDifferentials(Parameter parameter, double[] dRates, double[] dProportions);

    /**
     * @return true if the categories depend on the parameter through numerically computed quantiles, whose
     * truncation error (~1e-9 for the gamma shape) makes finite differences unreliable below steps of ~1e-4
     */
    default boolean hasQuantileNoise(Parameter parameter) {
        return false;
    }
}
