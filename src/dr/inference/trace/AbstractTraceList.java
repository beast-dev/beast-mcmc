/*
 * AbstractTraceList.java
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

package dr.inference.trace;

import java.util.List;

/**
 * @author Alexei Drummond
 */
public abstract class AbstractTraceList extends FilteredTraceList {
    /**
     * Please use {@link #getCorrelationStatistics(int) getCorrelationStatistics}
     * @param index
     * @return
     */
    @Deprecated
    public TraceDistribution getDistributionStatistics(int index) {
        return getCorrelationStatistics(index);
    }

    public TraceCorrelation getCorrelationStatistics(int index) {
        Trace trace = getTrace(index);
        if (trace != null)
            return trace.getTraceStatistics();
        return null;
    }

    public void analyseTrace(int index) {
        int start = (int) (getBurnIn() / getStepSize());

        Trace trace = getTrace(index);
        List<Double> values = trace.getValues(start, trace.getValueCount(), super.filtered);

        if (trace.getTraceType() == TraceType.CATEGORICAL) {
            trace.setTraceStatistics(new TraceCorrelation(values, trace.getCategoryLabelMap(), trace.getCategoryOrder(), getStepSize()));
        } else {
            trace.setTraceStatistics(new TraceCorrelation(values, trace.getTraceType(), getStepSize()));
        }

//        System.out.println("index = " + index + " :  " + trace.getName() + "     " + trace.getTraceType());
    }

//    abstract Trace getTrace(int index);

//    private TraceCorrelation[] traceStatistics = null;
}
