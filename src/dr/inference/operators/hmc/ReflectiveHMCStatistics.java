/*
 * ReflectiveHMCStatistics.java
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

package dr.inference.operators.hmc;

import dr.inference.loggers.LogColumn;
import dr.inference.loggers.Loggable;
import dr.inference.loggers.NumberColumn;

/**
 * Log columns for a reflective HMC operator. The counters are cumulative over the run, so per-sample values are the
 * differences between consecutive rows: wall events (time-0 events included), events that do not advance the drift,
 * walls (|J| summed over events), events with two or more walls, and rejections by cause (bounds check, zero-time cap,
 * maxEventsPerStep). largestWallSet and largestEventsInDrift are running maxima, not counters.
 */
public class ReflectiveHMCStatistics implements Loggable {

    private final ReflectiveHamiltonianMonteCarloOperator operator;
    private final String prefix;

    public ReflectiveHMCStatistics(ReflectiveHamiltonianMonteCarloOperator operator, String prefix) {
        this.operator = operator;
        this.prefix = prefix;
    }

    private interface Value {
        double get();
    }

    private LogColumn column(String name, final Value value) {
        return new NumberColumn(prefix + "." + name) {
            @Override
            public double getDoubleValue() {
                return value.get();
            }
        };
    }

    @Override
    public LogColumn[] getColumns() {
        return new LogColumn[]{
                column("events", operator::getEventCount),
                column("nonAdvancingEvents", operator::getNonAdvancingEventCount),
                column("walls", operator::getWallCount),
                column("multiWallEvents", operator::getMultiWallEventCount),
                column("largestWallSet", operator::getLargestWallSet),
                column("largestEventsInDrift", operator::getLargestEventsInDrift),
                column("boundsRejections", operator::getBoundsRejectionCount),
                column("zeroTimeCapRejections", operator::getZeroTimeCapCount),
                column("maxEventsRejections", operator::getMaxEventsRejectionCount),
        };
    }
}
