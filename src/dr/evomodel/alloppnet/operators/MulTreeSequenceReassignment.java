
/*
 * MulTreeSequenceReassignment.java
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

import dr.evomodel.alloppnet.speciation.MulSpeciesBindings;
import dr.evomodel.alloppnet.speciation.MulSpeciesTreeModel;
import dr.evomodel.alloppnet.parsers.MulTreeSequenceReassignmentParser;
import dr.inference.operators.SimpleMCMCOperator;
import dr.math.MathUtils;


/**
 * Changes sequence assignments for one gene, one or more individuals.
 * Very similar to AlloppSequenceReassignment.
 * 
 * @author Graham Jones
 *         Date: 20/12/2011
 */
public class MulTreeSequenceReassignment extends SimpleMCMCOperator {

	private final MulSpeciesTreeModel multree;
	private final MulSpeciesBindings mulspb;


	public MulTreeSequenceReassignment(MulSpeciesTreeModel multree, MulSpeciesBindings mulspb, double weight) {
		this.multree = multree;
		this.mulspb = mulspb;
		setWeight(weight);
	}	


	public String getPerformanceSuggestion() {
		return "None";
	}

	@Override
	public String getOperatorName() {
		return MulTreeSequenceReassignmentParser.MULTREE_SEQUENCE_REASSIGNMENT + "(" + multree.getId() +
		"," + mulspb.getId() + ")";
	}

	@Override
	public double doOperation() {
		multree.beginTreeEdit();
		if (MathUtils.nextInt(2) == 0) {
			mulspb.permuteOneSpeciesOneIndivForOneGene();
		} else {
			mulspb.permuteSetOfIndivsForOneGene();
		}
		
		multree.endTreeEdit();
		return 0;
	}

}


