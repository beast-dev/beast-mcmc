/*
 * JukesCantorDistanceMatrix.java
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

package dr.evolution.distance;

import dr.evolution.alignment.PatternList;

/**
 * compute jukes-cantor corrected distance matrix
 *
 * @version $Id: JukesCantorDistanceMatrix.java,v 1.4 2005/05/24 20:25:56 rambaut Exp $
 *
 * @author Andrew Rambaut
 * @author Korbinian Strimmer
 */
public class JukesCantorDistanceMatrix extends DistanceMatrix
{
	//
	// Public stuff
	//
	
	/** constructor */
	public JukesCantorDistanceMatrix()
	{
		super();
	}

	/** constructor taking a pattern source */
	public JukesCantorDistanceMatrix(PatternList patterns)
	{
		super(patterns);
	}

	/**
	 * set the pattern source
	 */
	public void setPatterns(PatternList patterns) {
		super.setPatterns(patterns);

        final int stateCount = patterns.getStateCount();

		const1 = ((double) stateCount - 1) / stateCount;
		const2 = ((double) stateCount) / (stateCount - 1) ;
	}
		
	/**
	 * Calculate a pairwise distance
	 */
	protected double calculatePairwiseDistance(int i, int j) {
		final double obsDist = super.calculatePairwiseDistance(i, j);
		
		if (obsDist == 0.0) return 0.0;
	
		if (obsDist >= const1) {
			return MAX_DISTANCE;
		} 
        
		final double expDist = -const1 * Math.log(1.0 - (const2 * obsDist));

		if (expDist < MAX_DISTANCE) {
			return expDist;
		} else {
			return MAX_DISTANCE;
		}
	}
	
	//
	// Private stuff
	//

    //used in correction formula
	private double const1, const2;
}
