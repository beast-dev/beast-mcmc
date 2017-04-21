/*
 * MoveLinkageGroup.java
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

import dr.evolution.util.Taxon;
import dr.evomodel.tree.HiddenLinkageModel;
import dr.evomodelxml.operators.MoveLinkageGroupParser;
import dr.inference.operators.SimpleMCMCOperator;
import dr.math.MathUtils;

/**
 * @author Aaron Darling
 */
public class MoveLinkageGroup extends SimpleMCMCOperator {

	HiddenLinkageModel hlm;	
	int readCount;
	int groupCount;
	
	public MoveLinkageGroup(HiddenLinkageModel hlm, double weight){
		this.hlm = hlm;
		readCount = hlm.getData().getReadsTaxa().getTaxonCount();
		groupCount = hlm.getLinkageGroupCount();
        setWeight(weight);
	}

	public double doOperation() {
		// pick a read uniformly at random, add it to a linkage group uniformly at random
		int r = MathUtils.nextInt(readCount);		
		Taxon read = hlm.getData().getReadsTaxa().getTaxon(r);
		int g=hlm.getLinkageGroupId(read);

		int newGroup = MathUtils.nextInt(groupCount);
		hlm.moveReadGroup(read, g, newGroup);

		// moves are symmetric -- same forward and backward proposal probabilities
		double logHastings = 0.0;
		return logHastings;
	}

	public String getOperatorName() {
		return MoveLinkageGroupParser.MOVE_LINKAGE_GROUP + "(" + hlm.getId() + ")";
	}

	public String getPerformanceSuggestion() {
		return "Ask Aaron Darling to write a better operator";
	}

}
