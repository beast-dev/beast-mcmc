/*
 * CombinedTraces.java
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

package dr.app.tracer.traces;

import dr.inference.trace.*;


/**
 * A class for analysing multiple tracesets
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: CombinedTraces.java,v 1.3 2006/11/29 16:04:12 rambaut Exp $
 */

public class CombinedTraces extends FilteredTraceList {

    public CombinedTraces(String name, TraceList[] traceLists) throws TraceException {

        if (traceLists == null || traceLists.length < 1) {
            throw new TraceException("Must have at least 1 Traces object in a CombinedTraces");
        }

        this.name = name;
        this.traceLists = new TraceList[traceLists.length];
        this.traceLists[0] = traceLists[0];

        for (int i = 1; i < traceLists.length; i++) {
            if (traceLists[i].getTraceCount() != traceLists[0].getTraceCount()) {
                throw new TraceException("Cannot add to a CombinedTraces: the count of new traces do not match existing traces");
            }

            if (traceLists[i].getStepSize() != traceLists[0].getStepSize()) {
                throw new TraceException("Cannot add to a CombinedTraces: the step sizes of the new traces do not match existing traces");
            }

            for (int j = 0; j < traceLists[0].getTraceCount(); j++) {
                if (!traceLists[i].getTraceName(j).equals(traceLists[0].getTraceName(j))) {
                    throw new TraceException("Cannot add to a CombinedTraces: new traces do not match existing traces");
                }
            }
            this.traceLists[i] = traceLists[i];
        }


    }

    /**
     * @return the name of this trace list
     */
    public String getName() {
        return name;
    }

    /**
     * @return the number of traces in this trace list
     */
    public int getTraceCount() {
        return traceLists[0].getTraceCount();
    }

    /**
     * @return the index of the trace with the given name
     */
    public int getTraceIndex(String name) {
        return traceLists[0].getTraceIndex(name);
    }

    /**
     * @return the name of the trace with the given index
     */
    public String getTraceName(int index) {
        return traceLists[0].getTraceName(index);
    }

    /**
     * @return the number of states excluding the burnin
     */
    public int getStateCount() {
        int sum = 0;
        for (TraceList traceList : traceLists) {
            sum += traceList.getStateCount();
        }
        return sum;
    }

    /**
     * @return the number of states in the burnin
     */
    public int getBurninStateCount() {
        return 0;
    }


    public boolean isIncomplete() {
        return false;
    }

    public int getBurnIn() {
        return 0;
    }

    /**
     * @return the size of the step between states
     */
    public int getStepSize() {
        return traceLists[0].getStepSize();
    }

    /**
     * @return the last state in the chain
     */
    public int getMaxState() {
        return getStateCount() * getStepSize();
    }

    public <T> void getValues(int index, T[] destination) {
        int offset = 0;
        for (TraceList traceList : traceLists) {
            traceList.getValues(index, destination, offset);
            offset += traceList.getStateCount();
        }
    }

    public <T> void getValues(int index, T[] destination, int offset) {
        for (TraceList traceList : traceLists) {
            traceList.getValues(index, destination, offset);
            offset += traceList.getStateCount();
        }
    }

    public <T> void getBurninValues(int index, T[] destination) {
        throw new UnsupportedOperationException("getBurninValues is not a valid operation on CombinedTracers");
    }

    /**
     * @return the trace distribution statistic object for the given index
     */
    public TraceDistribution getDistributionStatistics(int index) {
        return getCorrelationStatistics(index);
    }

    /**
     * @return the trace correlation statistic object for the given index
     */
    public TraceCorrelation getCorrelationStatistics(int index) {
        if (traceStatistics == null) {
            return null;
            // this can happen if the ESS has not been calculated yet.
//	    throw new RuntimeException("No ESS for combined traces? This is not supposed to happen.");
        }

        return traceStatistics[index];
    }

    public void analyseTrace(int index) {
        // no offset: burnin is handled inside each TraceList we own and invisible to us.

        if (traceStatistics == null) {
            traceStatistics = new TraceCorrelation[getTraceCount()];
        }

        Trace trace = getTrace(index);

        if (trace != null) {
            if (trace.getTraceType() == Double.class) {
                Double values[] = new Double[getStateCount()];

                getValues(index, values);
                traceStatistics[index] = new TraceCorrelation(values, getStepSize());

            } else if (trace.getTraceType() == Integer.class) {
                Integer values[] = new Integer[getStateCount()];

                getValues(index, values);
                traceStatistics[index] = new TraceCorrelation(values, getStepSize());

            } else if (trace.getTraceType() == String.class) {
                String values[] = new String[getStateCount()];

                getValues(index, values);
                traceStatistics[index] = new TraceCorrelation(values, getStepSize());

            } else {
                throw new RuntimeException("Trace type is not recognized: " + trace.getTraceType());
            }
        }

    }

    public Trace getTrace(int index) {
        for (TraceList traceList : traceLists) {
            if (traceList.getTrace(index).getTraceType() != traceLists[0].getTrace(index).getTraceType()) {
                return null; // trace type not comparable
            }
        }
        return traceLists[0].getTrace(index);
    }

    /**
     * @return the number of trace lists that make up this combined
     */
    public int getTraceListCount() {
        return traceLists.length;
    }

    /**
     * @param index the index of the trace
     * @return the ith TreeTrace
     */
    public TraceList getTraceList(int index) {
        return traceLists[index];
    }

    //************************************************************************
    // private methods
    //************************************************************************

    private TraceList[] traceLists = null;

//    private TraceCorrelation[] traceStatistics = null;

    private String name;
}