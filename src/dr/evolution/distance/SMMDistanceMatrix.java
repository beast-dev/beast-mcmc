/*
 * SMMDistanceMatrix.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
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
 *
 */

package dr.evolution.distance;

import dr.evolution.alignment.PatternList;

/**
 * @author Chieh-Hsi Wu
 * Date: 31/07/2009
 * Time: 4:37:38 PM
 * This class is used to calculate the distance between different microsatellite alleles
 */
public class SMMDistanceMatrix extends DistanceMatrix{

   /**
    * constructor taking a pattern source
    * 
    * @param patterns   a pattern of a microsatellite locus
    */
    public SMMDistanceMatrix(PatternList patterns) {
        super(patterns);
    }

    protected double calculatePairwiseDistance(int taxon1, int taxon2) {

        int[] pattern = patterns.getPattern(0);
        int state1 = pattern[taxon1];
        
        int state2 = pattern[taxon2];
        double distance = 0.0;

        if (!dataType.isAmbiguousState(state1) && !dataType.isAmbiguousState(state2))
            distance = Math.abs(state1 - state2);

        return distance;
    }
}
