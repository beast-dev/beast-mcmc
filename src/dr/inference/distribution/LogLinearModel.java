/*
 * LogLinearModel.java
 *
 * Copyright (c) 2002-2016 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.inference.distribution;

import dr.inference.model.Parameter;

/**
 * @author Marc A. Suchard
 */
@Deprecated // GLM stuff is now in inference.glm - this is here for backwards compatibility temporarily
public class LogLinearModel extends GeneralizedLinearModel {

    public LogLinearModel(Parameter dependentParam) {
        super(dependentParam);
    }

    @Override
    public double[] getXBeta() {
        double[] xBeta = super.getXBeta();
        for(int i=0; i<xBeta.length; i++)
            xBeta[i] = Math.exp(xBeta[i]);
        return xBeta;
    }

    protected double calculateLogLikelihood(double[] beta) {
        throw new RuntimeException("Not yet implemented.");
    }

    protected double calculateLogLikelihood() {
        throw new RuntimeException("Not yet implemented.");
    }

    protected boolean confirmIndependentParameters() {
        return false;
    }

    public boolean requiresScale() {
        return false;
    }
}
