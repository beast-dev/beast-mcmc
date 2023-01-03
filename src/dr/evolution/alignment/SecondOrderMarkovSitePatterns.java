/*
 * SecondOrderMarkovSitePatterns.java
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

package dr.evolution.alignment;

import dr.evolution.datatype.DataType;
import dr.evolution.datatype.PairedDataType;

/**
 * @author Xiang Ji
 * @author Jeff Thorne
 * @author Marc A. Suchard
 */
public class SecondOrderMarkovSitePatterns extends UncertainSiteList {

    private final PairedDataType dataType;
    private final SitePatterns sitePatterns;

    public SecondOrderMarkovSitePatterns(DataType dataType,
                                         SitePatterns sitePatterns) {
        super(dataType, sitePatterns);
        this.dataType = new PairedDataType(sitePatterns.getDataType());
        this.sitePatterns = sitePatterns;
        setUncertainSitePatterns();
    }

    @Override
    public int getSiteCount() {
        return sitePatterns.getSiteCount();
    }

    @Override
    public int getPatternState(int taxonIndex, int patternIndex) {
        throw new RuntimeException("Does not support tip states.");
    }

    private void setUncertainSitePatterns() {

        for (int siteIndex = 0; siteIndex < sitePatterns.getSiteCount(); siteIndex++) {

            double[][] uncertainPattern = new double[sitePatterns.getTaxonCount()][];

            for (int i = 0; i < sitePatterns.getTaxonCount(); i++) {
                double[] partials = new double[dataType.getStateCount()];

                final int observedState = sitePatterns.getSitePattern(siteIndex)[i];

                fillPartial(partials, observedState);

                uncertainPattern[i] = partials;
            }
            addPattern(uncertainPattern);
        }
    }

    private void fillPartial(double[] partials, int observedState) {
        for (int previousState = 0; previousState < sitePatterns.getStateCount(); previousState++) {
            partials[dataType.getState(previousState, observedState)] = 1.0;
        }
    }

}
