/*
 * DifferentiableDemographicModel.java
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

package dr.evomodel.coalescent.demographicmodel;

import dr.inference.model.Parameter;

/**
 * A demographic model whose demographic-function arguments can be differentiated with respect to its parameters.
 *
 * @author Xiang Ji
 * @author Marc A. Suchard
 */
public interface DifferentiableDemographicModel {

    /**
     * @return true if getArgumentDerivatives accepts the (scalar) parameter
     */
    boolean isDifferentiableWrtParameter(Parameter parameter);

    /**
     * @return d argument_n / d parameter for each argument n of the function returned by getDemographicFunction
     * (see DemographicFunction.getArgument), at the current parameter values and time offset
     */
    double[] getArgumentDerivatives(Parameter parameter);
}
