/*
 * OnePGammaDistributionModelParser.java
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

package dr.inferencexml.distribution;

import dr.inference.distribution.GammaDistributionModel;
import dr.inference.distribution.ParametricDistributionModel;
import dr.inference.model.Parameter;

/**
 */
public class OnePGammaDistributionModelParser extends DistributionModelParser {
    
    public String getParserName() {
        return GammaDistributionModel.ONE_P_GAMMA_DISTRIBUTION_MODEL;
    }

    ParametricDistributionModel parseDistributionModel(Parameter[] parameters, double offset) {
        return new GammaDistributionModel(parameters[0]);
    }

    public String[] getParameterNames() {
        return new String[]{SHAPE};
    }

    public String getParserDescription() {
        return "A model of a one parameter gamma distribution.";
    }

    public boolean allowOffset() {
        return false;
    }

    public Class getReturnType() {
        return GammaDistributionModel.class;
    }
}
