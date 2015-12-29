/*
 * JackknifePatterns.java
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

package dr.evolution.alignment;

import dr.math.MathUtils;

/**
 * Provides 50% jackknife replicate patterns 
 *
 * @version $Id: JackknifePatterns.java,v 1.8 2005/05/24 20:25:55 rambaut Exp $
 *
 * @author Andrew Rambaut
 */
public class JackknifePatterns extends ResamplePatterns
{

 	/**
	 * Constructor
	 */
	public JackknifePatterns() {
	}

 	/**
	 * Constructor
	 */
	public JackknifePatterns(SiteList patterns) {
		setPatterns(patterns);
	}

	/**
	 * Perform a jackknife resampling of the patterns
	 */
	public void resamplePatterns() {
	
		int siteCount = patterns.getSiteCount();
		int oldPatternCount = patterns.getPatternCount();
					
		int[] siteIndices = MathUtils.shuffled(siteCount);
		int n = siteCount / 2;
	
		patternIndices = new int[oldPatternCount];
		weights = new double[oldPatternCount];
				
		int pattern;
		
		patternCount = 0;
		
		for (int i = 0; i < n; i++) {

			pattern = patterns.getPatternIndex(siteIndices[i]);
			
			int j = 0;
			for (j = 0; j < patternCount; j++) {
				if (patternIndices[j] == pattern) {
					break;
				}
			} 
			
			if (j < patternCount) {
				weights[j] += 1.0;
			} else {			
				patternIndices[j] = pattern;
				weights[j] = 1.0;
				patternCount++;
			}
		}

	}

}
