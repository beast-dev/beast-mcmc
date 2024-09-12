/*
 * LogLinearModel.java
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

package dr.inference.distribution;

import dr.evomodel.substmodel.LogAdditiveCtmcRateProvider;
import dr.inference.model.Parameter;

import java.util.List;

/**
 * @author Marc A. Suchard
 */

public class LogLinearModel extends GeneralizedLinearModel implements LogAdditiveCtmcRateProvider.Integrated {

    public LogLinearModel(Parameter dependentParam) {
        super(dependentParam);
    }

    public double[] getSuperXBeta() {
        return super.getXBeta();
    }
    
    @Override
    public double[] getXBeta() {
        double[] xBeta = super.getXBeta();
        for(int i = 0; i < xBeta.length; i++) {
            xBeta[i] = Math.exp(xBeta[i]);
        }
        return xBeta;
    }

    public Parameter getLogRateParameter() {
        throw new RuntimeException("Not yet implemented.");
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

    @Override
    public LogLinearModel factory(List<Parameter> oldIndependentParameter, List<Parameter> newIndependentParameter)  {
        LogLinearModel newGLM = new LogLinearModel(dependentParam);
        for (int i = 0; i < numRandomEffects; i++) {
            newGLM.addRandomEffectsParameter(randomEffects.get(i));
        }
        for (int i = 0; i < numIndependentVariables; i++) {
            Parameter currentIndependentParameter = independentParam.get(i);
            final int index = oldIndependentParameter.indexOf(currentIndependentParameter);
            if (index != -1) {
                newGLM.addIndependentParameter(newIndependentParameter.get(index), designMatrix.get(i), indParamDelta.get(i));
            } else {
                newGLM.addIndependentParameter(currentIndependentParameter, designMatrix.get(i), indParamDelta.get(i));
            }
        }
        return newGLM;
    }
}
