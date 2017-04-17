/*
 * SimpleDistanceMatrix.java
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

package dr.app.realtime;

import dr.evolution.alignment.PatternList;
import dr.evolution.distance.DistanceMatrix;

/**
 * @author Guy Baele
 */
public class SimpleDistanceMatrix extends DistanceMatrix
{

	/** constructor */
	public SimpleDistanceMatrix()
	{
		super();
	}

	/** constructor taking a pattern source */
	public SimpleDistanceMatrix(PatternList patterns)
	{
		super(patterns);
	}

	/**
	 * set the pattern source
	 */
	public void setPatterns(PatternList patterns) {
		super.setPatterns(patterns);

        final int stateCount = patterns.getStateCount();
	}
		
	/**
	 * Calculate a pairwise distance
	 */
	protected double calculatePairwiseDistance(int taxon1, int taxon2) {

        int state1, state2;
        int n = patterns.getPatternCount();
        int gapState = patterns.getDataType().getGapState();
        int[] pattern;
        double weight, distance;
        double sumDistance = 0.0;
        double sumWeight = 0.0;

        for (int i = 0; i < n; i++) {
            pattern = patterns.getPattern(i);

            state1 = pattern[taxon1];
            state2 = pattern[taxon2];
            weight = patterns.getPatternWeight(i);

            //only count if neither state is a gap
            if (state1 != gapState && state2 != gapState) {
                if (!dataType.isAmbiguousState(state1) && !dataType.isAmbiguousState(state2) && state1 != state2) {
                    sumDistance += weight;
                }
                sumWeight += weight;
            }

        }

        distance = sumDistance / sumWeight;

		return distance;
	}

}
