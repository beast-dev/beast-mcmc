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

import java.util.ArrayList;
import java.util.List;

import dr.app.beauti.components.ComponentFactory;

/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @author Walter Xie
 * @version $Id$
 */
public class ClockModelOptions extends ModelOptions {

	// Instance variables
    private final BeautiOptions options;
   
    private FixRateType rateOptionClockModel = FixRateType.FIX_MEAN; 
    private double meanRelativeRate = 1.0;

    public ClockModelOptions(BeautiOptions options) {    	
    	this.options = options;
               
        initGlobalClockModelParaAndOpers();
    }
    
    private void initGlobalClockModelParaAndOpers() {
    	
        createParameter("allClockRates", "All the relative rates regarding clock models");
    	
        createOperator("deltaAllClockRates", "deltaAllClockRates",
        		"Delta exchange operator for all the relative rates regarding clock models", "allClockRates",      		 
        		OperatorType.DELTA_EXCHANGE, 0.75, rateWeights);
    	
    }

    /**
     * return a list of parameters that are required
     *
     * @param params the parameter list
     */
    public void selectParameters(List<Parameter> params) {    	    	
//    	if (options.clockModelOptions.getRateOptionClockModel() == FixRateType.FIX_MEAN) {
// TODO       	
//        }
    }

    /**
     * return a list of operators that are required
     *
     * @param ops the operator list
     */
    public void selectOperators(List<Operator> ops) {
    	if (options.clockModelOptions.getRateOptionClockModel() == FixRateType.FIX_MEAN) {
    		Operator deltaOperator = getOperator("deltaAllClockRates");
    		ops.add(deltaOperator);
    	}
    }

	
    
    /////////////////////////////////////////////////////////////
    public FixRateType getRateOptionClockModel() {
		return rateOptionClockModel;
	}

	public void setRateOptionClockModel(FixRateType rateOptionClockModel) {
		this.rateOptionClockModel = rateOptionClockModel;
	}

	public void setMeanRelativeRate(double meanRelativeRate) {
		this.meanRelativeRate = meanRelativeRate;
	}

	public double getMeanRelativeRate() {
		return meanRelativeRate;
	}

	@Override
	public String getPrefix() {
		// TODO Auto-generated method stub
		return null;
	}


}
