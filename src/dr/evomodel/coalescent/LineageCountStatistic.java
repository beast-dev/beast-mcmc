/*
 * LineageCountStatistic.java
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

package dr.evomodel.coalescent;

import dr.evolution.coalescent.IntervalList;
import dr.inference.model.Statistic;

/**
 * The number of lineages in each interval, of an interval list or of a likelihood that builds its
 * own intervals (OldAbstractCoalescentLikelihood).
 */
public class LineageCountStatistic extends Statistic.Abstract {

    private final IntervalList intervalList;
    private final OldAbstractCoalescentLikelihood acl;

    public LineageCountStatistic(IntervalList intervalList) {
        this.intervalList = intervalList;
        this.acl = null;
    }

    public LineageCountStatistic(OldAbstractCoalescentLikelihood acl) {
        this.intervalList = null;
        this.acl = acl;
    }

    public int getDimension() {
        return intervalList != null ? intervalList.getIntervalCount() : acl.getIntervalCount();
    }

    public double getStatisticValue(int dim) {
        return intervalList != null ? intervalList.getLineageCount(dim) : acl.getLineageCount(dim);
    }

}
