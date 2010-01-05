/*
 * TraceList.java
 *
 * Copyright (C) 2002-2006 Alexei Drummond and Andrew Rambaut
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
    int getBurnIn();

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
    int getStepSize();

    /**
     * @return the last state in the chain
     */
    int getMaxState();

    boolean isIncomplete();

    /**
     * get the values of trace with the given index (without burnin)
     *
     * @param index       the index of trace
     * @param destination the array to copy values into
     */
    void getValues(int index, double[] destination);

    /**
     * get the values of trace with the given index (without burnin)
     *
     * @param index       the index of trace
     * @param destination the array to copy values into
     * @param offset      the start position for copying into the destination array
     */
    void getValues(int index, double[] destination, int offset);

    /**
     * get the values of the burnin of the trace
     *
     * @param index       the index of trace
     * @param destination the array to copy values into
     */
    void getBurninValues(int index, double[] destination);

    /**
     * @param traceIndex the index of the trace
     * @return the trace distribution statistic object for the given index
     */
    TraceDistribution getDistributionStatistics(int traceIndex);

    /**
     * @param traceIndex the index of the trace
     * @return the trace correlation statistic object for the given index
     */
    TraceCorrelation getCorrelationStatistics(int traceIndex);

    void analyseTrace(int index);

}