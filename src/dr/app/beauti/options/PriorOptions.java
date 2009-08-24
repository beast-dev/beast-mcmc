/*
 * Copyright (C) 2002-2009 Alexei Drummond and Andrew Rambaut
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
 * BEAST is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package dr.app.beauti.options;

import java.util.List;

import dr.app.beauti.enumTypes.FixRateType;
import dr.math.MathUtils;


/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @author Walter Xie
 * @version $Id$
 */
public class PriorOptions extends ModelOptions {

	// Instance variables
    private final BeautiOptions options;
   

    public PriorOptions(BeautiOptions options) {    	
    	this.options = options;
               
        
    }
       
    
    /**
     * return a list of parameters that are required
     *
     * @param params the parameter list
     */
    public void selectParameters(List<Parameter> params) {    	    	
        
        double growthRateMaximum = 1E6;
        double birthRateMaximum = 1E6;
        double substitutionRateMaximum = 100;
        double logStdevMaximum = 10;
        double substitutionParameterMaximum = 100;
                
        double[] rootAndRate = options.clockModelOptions.calculateInitialRootHeightAndRate(options.dataPartitions);
        double avgInitialRootHeight = rootAndRate[0];
        double avgInitialRate = rootAndRate[1];
        
        if (options.clockModelOptions.getRateOptionClockModel() == FixRateType.FIX_MEAN
    			|| options.clockModelOptions.getRateOptionClockModel() == FixRateType.RElATIVE_TO) {
        	
            growthRateMaximum = 1E6 * avgInitialRate;
            birthRateMaximum = 1E6 * avgInitialRate;   
        }
        
//        if (options.clockModelOptions.getRateOptionClockModel() == FixRateType.FIX_MEAN) {
//            double rate = options.clockModelOptions.getMeanRelativeRate();
//
//            growthRateMaximum = 1E6 * rate;
//            birthRateMaximum = 1E6 * rate;
//
//            if (options.hasData()) {
//                initialRootHeight = meanDistance / rate;
//
//                initialRootHeight = round(initialRootHeight, 2);
//            }
//
//        } else {
//            if (options.maximumTipHeight > 0) {
//                initialRootHeight = options.maximumTipHeight * 10.0;
//            }
//
//            initialRate = round((meanDistance * 0.2) / initialRootHeight, 2);
//        }

        double timeScaleMaximum = MathUtils.round(avgInitialRootHeight * 1000.0, 2);

        for (Parameter param : params) {
//            if (dataReset) param.priorEdited = false;

            if (!param.priorEdited) {
                switch (param.scale) {
                    case TIME_SCALE:
                        param.uniformLower = Math.max(0.0, param.lower);
                        param.uniformUpper = Math.min(timeScaleMaximum, param.upper);
                        param.initial = avgInitialRootHeight;
                        break;
                    case T50_SCALE:
                        param.uniformLower = Math.max(0.0, param.lower);
                        param.uniformUpper = Math.min(timeScaleMaximum, param.upper);
                        param.initial = avgInitialRootHeight / 5.0;
                        break;
                    case GROWTH_RATE_SCALE:
                        param.uniformLower = Math.max(-growthRateMaximum, param.lower);
                        param.uniformUpper = Math.min(growthRateMaximum, param.upper);
                        break;
                    case BIRTH_RATE_SCALE:
                        param.uniformLower = Math.max(0.0, param.lower);
                        param.uniformUpper = Math.min(birthRateMaximum, param.upper);
                        break;
                    case SUBSTITUTION_RATE_SCALE:
                        param.uniformLower = Math.max(0.0, param.lower);
                        param.uniformUpper = Math.min(substitutionRateMaximum, param.upper);
                        param.initial = avgInitialRate;
                        break;
                    case LOG_STDEV_SCALE:
                        param.uniformLower = Math.max(0.0, param.lower);
                        param.uniformUpper = Math.min(logStdevMaximum, param.upper);
                        break;
                    case SUBSTITUTION_PARAMETER_SCALE:
                        param.uniformLower = Math.max(0.0, param.lower);
                        param.uniformUpper = Math.min(substitutionParameterMaximum, param.upper);
                        break;

                    case UNITY_SCALE:
                        param.uniformLower = 0.0;
                        param.uniformUpper = 1.0;
                        break;

                    case ROOT_RATE_SCALE:
                        param.initial = avgInitialRate;
                        param.gammaAlpha = 0.5;
                        param.gammaBeta = param.initial / 0.5;
                        break;

                    case LOG_VAR_SCALE:
                        param.initial = avgInitialRate;
                        param.gammaAlpha = 2.0;
                        param.gammaBeta = param.initial / 2.0;
                        break;

                }
                if (param.isNodeHeight) { //TODO affecting "treeModel.rootHeight", need to review
                    param.lower = options.maximumTipHeight;
                    param.uniformLower = options.maximumTipHeight;
                    param.uniformUpper = timeScaleMaximum;
                    param.initial = avgInitialRootHeight;
                }
            }
        }

//        dataReset = false;
        
        
    }

    /**
     * return a list of operators that are required
     *
     * @param ops the operator list
     */
    public void selectOperators(List<Operator> ops) {
                
    }
    
    /////////////////////////////////////////////////////////////



}
