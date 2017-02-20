/*
 * MicrosatelliteModelSelectOperator.java
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

package dr.evomodel.operators;

import dr.inference.model.Parameter;
import dr.inference.operators.SimpleMCMCOperator;
import dr.math.MathUtils;

/**
 * @author Chieh-Hsi Wu
 *
 * Operator that selects a microsatellite model from a group provided by the user.
 */
public class MicrosatelliteModelSelectOperator extends SimpleMCMCOperator {
    private Parameter parameter;
    private Parameter[] indicators;
    public MicrosatelliteModelSelectOperator(Parameter parameter, Parameter[] indicators, double weight){
        this.parameter = parameter;
        this.indicators = indicators;
        setWeight(weight);
    }

    public String getOperatorName(){
        return "msatModelSelectOperator("+parameter.getParameterName()+")";
    }

    public final String getPerformanceSuggestion() {
        return "no suggestions available";
    }

    public double doOperation() {
        int index = MathUtils.nextInt(indicators.length);
        //System.out.println(index);
        Parameter newModel = indicators[index];
        for(int i = 0; i < parameter.getDimension() -1 ; i++){
            parameter.setParameterValueQuietly(i,newModel.getParameterValue(i));
        }
        parameter.setParameterValueNotifyChangedAll(
                parameter.getDimension()-1,
                newModel.getParameterValue(parameter.getDimension()-1)
        );
        //System.out.println(parameter+"; "+indicators[index]);
        return 0.0;
    }



}
