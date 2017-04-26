/*
 * MulSpeciesTreePrior.java
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

/**
 *
 * @author Graham Jones
 *         Date: 21/12/2011
 */

package dr.evomodel.alloppnet.speciation;

import dr.evomodel.speciation.SpeciationModel;
import dr.inference.model.Likelihood;



public class MulSpeciesTreePrior  extends Likelihood.Abstract {
	MulSpeciesTreeModel mulsptree;
	SpeciationModel prior;

	public MulSpeciesTreePrior(SpeciationModel prior, MulSpeciesTreeModel mulsptree) {
		super(prior);
		this.mulsptree = mulsptree;
		this.prior = prior;
	}

	
	@Override
	protected boolean getLikelihoodKnown() {
		return false;
	}

	
	
	@Override
	protected double calculateLogLikelihood() {
		double lhood = prior.calculateTreeLogLikelihood(mulsptree);
		return lhood;
	}
}
