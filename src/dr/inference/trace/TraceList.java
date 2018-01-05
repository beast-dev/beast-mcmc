/*
 * TraceList.java
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
 * An interface and default class that stores a set of traces from a single chain
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: TraceList.java,v 1.2 2006/05/22 15:52:06 rambaut Exp $
 */

public interface TraceList {

    /**
     * @return the name of this trace list
     */
    String getName();

    /**
     * @return the the full name of this trace list (possibly the path)
     */
    String getFullName();

    /**
     * @return the number of traces in this trace list
     */
    int getTraceCount();

    /**
     * @param name trace name
     * @return the index of the trace with the given name
     */
    int getTraceIndex(String name);

    /**
     * @param index trace ordinal index
     * @return the name of the trace with the given index
     */
    String getTraceName(int index);

    /**
     * @return the burn-in for this trace list (the number of actual states to discard)
     */
    long getBurnIn();

    /**
     * @return the number of states in the traces (without the burnin)
     */
    int getStateCount();

    /**
     * @return the number of states in the burnin
     */
    int getBurninStateCount();

    /**
     * @return the size of the step between states
     */
    long getStepSize();

    /**
     * @return the last state in the chain
     */
    long getMaxState();

    boolean isIncomplete();

    /**
     * get the values of trace with the given index (without burnin)
     * @param index       the index of trace
     * @param fromIndex   low endpoint (inclusive) of the subList.
     * @param toIndex     high endpoint (exclusive) of the subList.
     * @return The list of values (which are selected values if filter applied)
     */
    List<Double> getValues(int index, int fromIndex, int toIndex);

    /**
     * get the values of trace with the given index (without burnin)
     * @param index       the index of trace
     * @return The list of values (which are selected values if filter applied)
     */
    List<Double> getValues(int index);

    /**
     * get the values of the burnin of the trace
     * @param index       the index of trace
     * @return The list of values (which are selected values if filter applied)
     */
    List<Double> getBurninValues(int index);

    /**
     * @param traceIndex the index of the trace
     * @return the trace distribution statistic object for the given index
     */
    @Deprecated
    TraceDistribution getDistributionStatistics(int traceIndex);

    /**
     * @param traceIndex the index of the trace
     * @return the trace correlation statistic object for the given index
     */
    TraceCorrelation getCorrelationStatistics(int traceIndex);

    // create TraceCorrelation regarding Trace
    void analyseTrace(int index);

    Trace getTrace(int index);

//    public interface D extends TraceList {
//        Double[] getValues(int index, int length);
//        Double[] getValues(int index, int length, int offset);
//        Double[] getBurninValues(int index, int length);
//    }
}