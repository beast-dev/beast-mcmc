/*
 * BootstrapPatterns.java
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

package dr.evolution.alignment;

import dr.math.MathUtils;

/**
 * Provides bootstrap replicate patterns
 *
 *
 * @author Andrew Rambaut
 */
public class BootstrapPatterns extends ResamplePatterns {

 	/**
	 * Constructor
	 */
	public BootstrapPatterns() {
	}

 	/**
	 * Constructor
	 */
	public BootstrapPatterns(SiteList patterns) {
		setPatterns(patterns);
	}

	/**
	 * Perform a bootstrap resampling of the patterns
	 */
	public void resamplePatterns() {

		int siteCount = patterns.getSiteCount();
		int oldPatternCount = patterns.getPatternCount();

		patternIndices = new int[oldPatternCount];
		weights = new double[oldPatternCount];

		int site, pattern;

		patternCount = 0;

		for (int i = 0; i < siteCount; i++) {
			site = MathUtils.nextInt(siteCount);
			pattern = patterns.getPatternIndex(site);

			int j = 0;
			for (j = 0; j < patternCount; j++) {
				if (patternIndices[j] == pattern) {
					break;
				}
			}

			if (j < patternCount) {
				weights[j] += 1.0;
			} else {
				patternIndices[patternCount] = pattern;
				weights[patternCount] = 1.0;
				patternCount++;
			}
		}

	}
}
