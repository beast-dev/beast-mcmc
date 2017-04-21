/*
 * AlloppNetworkPriorModel.java
 *
 * Copyright (c) 2002-2017 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodel.alloppnet.speciation;

import dr.evolution.util.Units;
import dr.evomodel.alloppnet.parsers.AlloppNetworkPriorModelParser;
import dr.inference.distribution.ParametricDistributionModel;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.inference.model.Variable.ChangeType;



/**
 * 
 * Model for prior of an allopolyploid network.
 * 
 * @author Graham Jones
 *         Date: 01/07/2011
 */


/*





 */

public class AlloppNetworkPriorModel extends AbstractModel  implements Units {
    private Parameter rate;
    private Parameter popsf;
    private ParametricDistributionModel tpmodel;
    private ParametricDistributionModel rpmodel;
    private ParametricDistributionModel hpmodel;
    private Units.Type units;


    public AlloppNetworkPriorModel(Parameter eventrate, Parameter popscalingfactor,
                            ParametricDistributionModel tippopmodel, ParametricDistributionModel rootpopmodel,
                            ParametricDistributionModel hybpopmodel, Units.Type  units) {
        super(AlloppNetworkPriorModelParser.ALLOPPNETWORKPRIORMODEL);
        this.rate = eventrate;
        this.popsf = popscalingfactor;
        this.tpmodel = tippopmodel;
        this.rpmodel = rootpopmodel;
        this.hpmodel = hybpopmodel;
        this.units = units;
        addVariable(rate);
        addVariable(popsf);
    }


    Parameter getRate() {
		return rate;
	}


    Parameter getPopScalingFactor() {
        return popsf;
    }

    ParametricDistributionModel getTipPopModel() {
        return tpmodel;
    }

    ParametricDistributionModel getRootPopModel() {
        return rpmodel;
    }

    ParametricDistributionModel getHybridPopModel() {
        return hpmodel;
    }


    @Override
	protected void handleModelChangedEvent(Model model, Object object, int index) {
		if (AlloppSpeciesNetworkModel.DBUGTUNE)
			System.out.println("AlloppNetworkPriorModel.handleModelChangedEvent() " + model.getId());
	}

	@Override
	protected void handleVariableChangedEvent(Variable variable, int index,
			ChangeType type) {
		if (AlloppSpeciesNetworkModel.DBUGTUNE)
			System.out.println("AlloppNetworkPriorModel.handleModelChangedEvent() " + variable.getId());
	}

	@Override
	protected void storeState() {
		// addVariable(rate), addVariable(popsf) deal with this
		if (AlloppSpeciesNetworkModel.DBUGTUNE)
			System.out.println("AlloppNetworkPriorModel.storeState()");

	}

	@Override
	protected void restoreState() {
		// addVariable(rate), addVariable(popsf) deal with this
		if (AlloppSpeciesNetworkModel.DBUGTUNE)
			System.out.println("AlloppNetworkPriorModel.restoreState()");

	}

	@Override
	protected void acceptState() {
	}
	
	public Type getUnits() {
		return units;
	}

	public void setUnits(Type units) {
		// grjtodo-oneday allow units other than substitutions
		
	}	

}
