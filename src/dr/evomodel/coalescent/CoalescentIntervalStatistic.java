/*
 * CoalescentIntervalStatistic.java
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

package dr.evomodel.coalescent;

import dr.inference.model.Statistic;
import dr.math.Binomial;


public class CoalescentIntervalStatistic extends Statistic.Abstract {

    private final CoalescentIntervalProvider coalescent;
    private final boolean rescaleToNe;

    public CoalescentIntervalStatistic(CoalescentIntervalProvider coalescent) {
        this(coalescent, false);
    }

    public CoalescentIntervalStatistic(CoalescentIntervalProvider coalescent, boolean rescaleToNe) {
        this.coalescent = coalescent;
        this.rescaleToNe = rescaleToNe;
    }

    public int getDimension() {

        throw new RuntimeException("the use of CoalescentIntervalStatistic has been deprecated");

        //return coalescent.getCoalescentIntervalDimension();

    }

    public double getStatisticValue(int i) {

        throw new RuntimeException("the use of CoalescentIntervalStatistic has been deprecated");

        /*double interval = coalescent.getCoalescentInterval(i);

        if (rescaleToNe) {
            int lineages = coalescent.getCoalescentIntervalLineageCount(i);
            interval *= Binomial.choose2(lineages);
            // TODO Double-check; maybe need to return 1/interval or divide by choose2(lineages)
        }*/

        //return interval;

    }
    
    public String getStatisticName() {
    	return "coalescentIntervalStatistic";
    }
    
}
