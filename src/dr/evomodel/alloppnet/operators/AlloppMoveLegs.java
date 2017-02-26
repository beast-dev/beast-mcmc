/*
 * AlloppMoveLegs.java
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

package dr.evomodel.alloppnet.operators;

import dr.evomodel.alloppnet.speciation.AlloppSpeciesBindings;
import dr.evomodel.alloppnet.speciation.AlloppSpeciesNetworkModel;
import dr.evomodel.alloppnet.parsers.AlloppMoveLegsParser;
import dr.inference.operators.SimpleMCMCOperator;


/**
 * 
 * @author Graham Jones
 *         Date: 31/08/2011
 */


// 2012-07 ood
public class AlloppMoveLegs extends SimpleMCMCOperator {

	private final AlloppSpeciesNetworkModel apspnet;
	private final AlloppSpeciesBindings apsp;
	
	public AlloppMoveLegs(AlloppSpeciesNetworkModel apspnet, AlloppSpeciesBindings apsp, double weight) {
		this.apspnet = apspnet;
		this.apsp = apsp;
		setWeight(weight);
	}	

	
	public String getPerformanceSuggestion() {
		return "None";
	}

	@Override
	public String getOperatorName() {
		return AlloppMoveLegsParser.MOVE_LEGS + "(" + apspnet.getId() +
		"," + apsp.getId() + ")";
	}

	@Override
	public double doOperation() {
		apspnet.beginNetworkEdit();
		apspnet.moveLegs();		
		apspnet.endNetworkEdit();
		return 0;
	}

}
