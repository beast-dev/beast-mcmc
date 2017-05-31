/*
 * ShuffledSiteList.java
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
 * Provides bootstrap replicate patterns 
 *
 * @version $Id: ShuffledSiteList.java,v 1.3 2005/04/20 21:26:19 rambaut Exp $
 *
 * @author Andrew Rambaut
 */
public class ShuffledSiteList extends ResamplePatterns implements SiteList
{

 	/**
	 * Constructor
	 */
	public ShuffledSiteList() {
	}

 	/**
	 * Constructor
	 */
	public ShuffledSiteList(SiteList patterns) {
		setPatterns(patterns);
	}

    public void setPatterns(SiteList patterns) {
        this.patterns = patterns;
        resamplePatterns();
    }
	/**
	 * Perform a bootstrap resampling of the patterns
	 */
	public void resamplePatterns() {
	
		int siteCount = patterns.getSiteCount();
        siteIndices = MathUtils.shuffled(siteCount);
	}

    public int getSiteCount() {
        return siteIndices.length;
    }

    public int[] getSitePattern(int siteIndex) {
        return patterns.getSitePattern(siteIndices[siteIndex]);
    }

    @Override
    public double[][] getUncertainSitePattern(int siteIndex) {
        return patterns.getUncertainSitePattern(siteIndices[siteIndex]);
	}

    public int getPatternIndex(int siteIndex) {
        return patterns.getPatternIndex(siteIndices[siteIndex]);
    }

    public int getState(int taxonIndex, int siteIndex) {
        return patterns.getState(taxonIndex, siteIndices[siteIndex]);
    }

    @Override
    public double[] getUncertainState(int taxonIndex, int siteIndex) {
        return patterns.getUncertainState(taxonIndex, siteIndices[siteIndex]);
	}

    int siteIndices[] = null;

    @Override
    public boolean areUncertain() {
        return false;
    }
}
