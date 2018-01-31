/*
 * LogFileTraces.java
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

import java.io.*;
import java.util.ArrayList;
import java.util.List;
import java.util.StringTokenizer;

/**
 * A class that stores a set of traces from a single chain
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: LogFileTraces.java,v 1.4 2006/11/30 17:39:29 rambaut Exp $
 */

public class LogFileTraces extends AbstractTraceList {

    public LogFileTraces(String name, File file) {
        this.name = name;
        this.file = file;
    }

    /**
     * @return the name of this traceset
     */
    public String getName() {
        return name;
    }

    /**
     * @return the path of this traceset
     */
    public String getFullName() {
        return file.getPath();
    }

    public File getFile() {
        return file;
    }

    /**
     * @return the last state in the chain
     */
    public long getMaxState() {
        return lastState;
    }

    public boolean isIncomplete() {
        return false;
    }

    /**
     * @return the number of states excluding the burnin
     */
    public int getStateCount() {
        // This is done as two integer divisions to ensure the same rounding for
        // the burnin...
        return (int) (((lastState - firstState) / stepSize) - (getBurnIn() / stepSize) + 1);
    }

    /**
     * @return the number of states in the burnin
     */
    public int getBurninStateCount() {
        return (int) (getBurnIn() / stepSize);
    }

    /**
     * @return the size of the step between states
     */
    public long getStepSize() {
        return stepSize;
    }

    public long getBurnIn() {
        return burnIn;
    }

    /**
     * @return the number of traces in this traceset
     */
    public int getTraceCount() {
        return traces.size();
    }

    /**
     * @return the index of the trace with the given name
     */
    public int getTraceIndex(String name) {
        for (int i = 0; i < traces.size(); i++) {
            Trace trace = getTrace(i);
            if (name.equals(trace.getName())) {
                return i;
            }
        }
        return -1;
    }

    /**
     * @return the name of the trace with the given index
     */
    public String getTraceName(int index) {
        return getTrace(index).getName();
    }

    /**
     * @param index requested trace index
     * @return the trace for a given index
     */
    public Trace getTrace(int index) {
        return traces.get(index);
    }

    public void setBurnIn(long burnin) {
        this.burnIn = burnin;
        for (Trace trace : traces) {
            trace.setTraceStatistics(null);
        }
    }

    public double getStateValue(int trace, int index) {
        return (Double) getTrace(trace).getValue(index + (int) (getBurnIn() / stepSize));
    }

    /**
     * Read several consecutive values of one state into a destination array
     *
     * @param nState      State index number
     * @param destination array to store result
     * @param offset      first trace index
     */
    public void getStateValues(int nState, double[] destination, int offset) {
        final int index1 = nState + (int) (getBurnIn() / stepSize);
        for (int k = 0; k < destination.length; ++k) {
            destination[k] = (Double) getTrace(k + offset).getValue(index1);
        }
    }

    /**
     * Use the flag boolean[] filtered in FilteredTraceList
     * to determine whether to remove filtered values,
     * when filtered != null
     * @param index       the index of trace
     * @param fromIndex   low endpoint (inclusive) of the subList.
     * @param toIndex     high endpoint (exclusive) of the subList.
     * @return
     */
    public List<Double> getValues(int index, int fromIndex, int toIndex) {
        List<Double> values = null;
        try {
            Trace trace = getTrace(index);
            values = trace.getValues(fromIndex, toIndex, super.filtered);
        } catch (IndexOutOfBoundsException e) {
            System.err.println("getValues error: trace index = " + index);
        }
        return values;
    }

    public List<Double> getValues(int index) {
        return this.getValues(index, getBurninStateCount(), getTrace(index).getValueCount());
    }

    public List<Double> getBurninValues(int index) {
        return this.getValues(index, 0, getBurninStateCount());
    }

    /**
     * Use the {@link #loadTraces(File) loadTraces} method,
     * where <code>File</code> is defined from the constructor.
     *
     * @throws TraceException
     * @throws IOException
     */
    public void loadTraces() throws TraceException, IOException {
        loadTraces(file);
    }

    /**
     * Read through <code>File</code> created from a log file,
     * fill in <code>traces</code> list, and set <code>TraceType</code>.
     *
     * @param file <code>File</code>
     * @throws TraceException
     * @throws IOException
     */
    public void loadTraces(File file) throws TraceException, IOException {
        final Reader reader = new FileReader(file);
        loadTraces(reader);
        reader.close();
    }

    /**
     * Read through <code>InputStream</code> created from a log file,
     * fill in <code>traces</code> list, and set <code>TraceType</code>.
     *
     * @param in <code>InputStream</code>
     * @throws TraceException
     * @throws IOException
     */
    public void loadTraces(InputStream in) throws TraceException, IOException {
        final Reader reader = new InputStreamReader(in);
        loadTraces(reader);
        reader.close();
    }

    /**
     * Read through either <code>FileReader</code> or <code>InputStreamReader</code>
     * created from a log file,
     * fill in <code>traces</code> list, and set <code>TraceType</code>.
     *
     * @param r The input for <code>TrimLineReader</code>.
     *          Use either <code>FileReader</code> or <code>InputStreamReader</code>
     * @throws TraceException
     * @throws java.io.IOException
     */
    private void loadTraces(Reader r) throws TraceException, java.io.IOException {

        final TrimLineReader reader = new LogFileTraces.TrimLineReader(r);

        // Read through to first token
        StringTokenizer tokens = reader.tokenizeLine();

        if (tokens == null) {
            throw new TraceException("Trace file is empty.");
        }

        // read over empty lines
        while (!tokens.hasMoreTokens()) {
            tokens = reader.tokenizeLine();
        }

        // skip the first column which should be the state number
        String token = tokens.nextToken();

        // lines starting with [ are ignored, assuming comments in MrBayes file
        // lines starting with # are ignored, assuming comments in Migrate or BEAST file
        while (token.startsWith("[") || token.startsWith("#")) {
            // readTraceType(token, tokens); // using # to define type
            tokens = reader.tokenizeLine();

            // read over empty lines
            while (!tokens.hasMoreTokens()) {
                tokens = reader.tokenizeLine();
            }

            // read state token and ignore
            token = tokens.nextToken();
        }

        // read label tokens
        String[] labels = new String[tokens.countTokens()];

        for (int i = 0; i < labels.length; i++) {
            labels[i] = tokens.nextToken();
            addTrace(labels[i]);
        }

        int traceCount = getTraceCount();

        long num_samples = 0;

        String line = reader.readLine();
        tokens = reader.getStringTokenizer(line);
        String lastLine = line;
        while (tokens != null && tokens.hasMoreTokens()) {

            String stateString = tokens.nextToken();
            long state = 0;

            try {
                try {
                    // Changed this to parseDouble because LAMARC uses scientific notation for the state number
                    state = (long) Double.parseDouble(stateString);
                } catch (NumberFormatException nfe) {
                    throw new TraceException("Unable to parse state number in column 1 (Line " +
                            reader.getLineNumber() + ")");
                }

                if (num_samples < 1) {
                    // MrBayes puts 1 as the first state, BEAST puts 0
                    // In order to get the same gap between subsequent samples,
                    // we force this to 0.
                    if (state == 1) state = 0;
                }
                num_samples += 1;

                if (!addState(state, num_samples)) {
                    throw new TraceException("State " + state + " is not consistent with previous spacing (Line " +
                            reader.getLineNumber() + ")");
                }

            } catch (NumberFormatException nfe) {
                throw new TraceException("State " + state + ":Expected real value in column " + reader.getLineNumber());
            }

            for (int i = 0; i < traceCount; i++) {
                if (tokens.hasMoreTokens()) {
                    String value = tokens.nextToken();

                    addParsedValue(i, value);

                } else {
                    throw new TraceException("State " + state + ": missing values at line " + reader.getLineNumber());
                }
            }
            // used to keep the last valid line
            lastLine = line;
//            tokens = reader.tokenizeLine();
            line = reader.readLine();
            tokens = reader.getStringTokenizer(line);
        }

        if (num_samples == 0)
            throw new TraceException("Incorrect file format, no sample is found !");

        burnIn =  lastState / 10;

        if (lastState < 0)
            lastState = firstState;
        if (stepSize < 0 && lastState > 0)
            stepSize = lastState;

        validateTraceType(lastLine);
        validateUniqueValues();
    }
    public static final int MIN_SAMPLE = 5; // used in StatisticsModel

    private final int MAX_UNIQUE_VALUE = 200;
    // change integer type into real, if too many unique values
    private void validateUniqueValues() throws TraceException {
        for (int id = 0; id < getTraceCount(); id++) {
            Trace trace = getTrace(id);
            if (trace.getTraceType().isInteger()) {
                int uniqueValue = trace.getUniqueValueCount();
                if (uniqueValue > MAX_UNIQUE_VALUE) {
                    System.out.println("Too many unique values (>" + MAX_UNIQUE_VALUE +
                            ") found in trace " + trace.getName() + " at " + id);
                    changeTraceType(id, TraceType.REAL);
                }
            }
        }
    }


    // validate TraceType at the last value of trace,
    // in case integer is logged for double values in the first (even several) row.
    // it must use original line, because the data type of values in traces are changed
    private void validateTraceType(String lastLine) throws TraceException {
        String[] values = lastLine.split("\\t");
        // the 1st is state
        for (int i=1; i < values.length; i++) {
            int traceId = i-1;
            Trace trace = getTrace(traceId);
            // avoid to assign integer to double incorrectly
            if (trace.getTraceType().isInteger() && NumberUtils.hasDecimalPoint(values[i]))
                changeTraceType(traceId, TraceType.REAL);
        }
    }

    /**
     * Auto assign INTEGER or CATEGORICAL type to traces
     * according their values in the first line.
     * Default type is REAL.
     *
     * @param traceIndex
     * @param value
     */
    private void addParsedValue(int traceIndex, String value) throws TraceException {
        Trace trace = getTrace(traceIndex);

        TraceType type = trace.getTraceType();

        if (type != TraceType.CATEGORICAL) {
            // once categorical, always categorical

            if (NumberUtils.isNumber(value)) {
                if (type != TraceType.REAL) {
                    // once real, always real (unless categorical)

                    if (NumberUtils.hasDecimalPoint(value)) { // Real
                        // if a single number has a decimal point then switch to real
                        type = TraceType.REAL;
//                        System.out.println("Auto detect " + type + " type for trace " + name + " at " + traceIndex);
                        changeTraceType(traceIndex, type);
                    }
                }
            } else { // String
                type = TraceType.CATEGORICAL;
//                System.out.println("Auto detect " + type + " type for trace " + name + " at " + traceIndex);
                trace = changeTraceType(traceIndex, type);
            }
        }

        if (type == TraceType.REAL) {
            trace.add(Double.parseDouble(value));
        } else if (type == TraceType.INTEGER) {
            trace.add(Double.parseDouble(value));
        } else if (type == TraceType.CATEGORICAL) {
            trace.add(value);
        } else {
            throw new IllegalArgumentException("Unsupported TraceType");
        }

    }

    //************************************************************************
    // private methods
    //************************************************************************

    // These methods are used by the load function, above

    /**
     * Add a trace for a statistic of the given name
     *
     * @param name trace name
     */
    private void addTrace(String name) {
        // Start by assuming integer. This will be switched if added data
        // contradicts.
        traces.add(new Trace(name, TraceType.INTEGER));
    }

    /**
     * Add a trace for a statistic of the given name
     *
     * @param name trace name
     */
    private void addTrace(String name, TraceType traceType) {
        traces.add(new Trace(name, traceType));
    }

    public Trace changeTraceType(int id, TraceType newType) throws TraceException {
        if (id >= getTraceCount() || id < 0)
            throw new TraceException("Invalid trace id : " + id + ", which should 0 < and >= " + getTraceCount());
        Trace trace = getTrace(id);
        TraceType oldType = trace.getTraceType();
        if (oldType != newType) {
            if (oldType.isCategorical()) {
                throw new TraceException("A categorical type cannot be changed to anything else.");
            }

            Trace newTrace = new Trace(trace.getName(), newType);

            if (newType.isDiscrete()) {
                int uniqueValue = trace.getUniqueValueCount();
                if (uniqueValue > MAX_UNIQUE_VALUE)
                    throw new TraceException("Type change is failed, because too many unique values (>" +
                            MAX_UNIQUE_VALUE + ") are found !");
            }

            if (newType.isCategorical()) {
                try {
                    // include integer trace
                    for (int i = 0; i < trace.getValueCount(); i++) {
                        // use add(String value) to set categoryLabelMap ...
                        newTrace.add(Double.toString(trace.getValue(i)));
                    }
                } catch (Exception e) {
                    throw new TraceException("Type change is failed, when parsing " + oldType +
                            " to " + newType + " in trace " + trace.getName());
                }

                if (newTrace.getValueCount() != trace.getValueCount())
                    throw new TraceException("Type change is failed, because values size is different after copy !");

                traces.set(id, newTrace);

            } else {
                trace.setTraceType(newType); // change between numeric
            }

            // do not need code below, they are set by 'add(String value)',
            // and they will overwrite the correct values
//            // copy the categorical values across in case it is switched back
//            newTrace.categoryValueList = trace.categoryValueList;
//            newTrace.categoryLabelMap = trace.categoryLabelMap;
//            newTrace.orderType = trace.orderType;
//            newTrace.categoryOrder = trace.categoryOrder;
//            newTrace.uniqueValues = trace.uniqueValues;

//            System.out.println("Change " + oldType + " to " + newType + " type for trace " + trace.getName() + " at " + id);
            return newTrace;
        }
        return trace;
    }

    /**
     * Add a state number for these traces. This should be
     * called before adding values for each trace. The spacing
     * between stateNumbers should remain constant.
     *
     * @param stateNumber the state
     * @param num_samples the number of samples (rows)
     * @return false if the state number is inconsistent
     */
    private boolean addState(long stateNumber, long num_samples) {
        if (firstState < 0) { // it can use num_samples==1 to replace firstState < 0
            firstState = stateNumber;
        } else if (secondState < 0) {
            secondState = stateNumber;
        } else if (stepSize < 0) {
            // delay setting of the stepSize until the step between
            // the second and third step in case the first step is
            // 1 (i.e., MrBayes) and the stepsize is 1.
            stepSize = stateNumber - secondState;
        } else {
            long step = stateNumber - lastState;
            if (step != stepSize) {
                System.out.println("step: " + step + " != " + stepSize);
                return false;
            }
        }
//        System.out.println(num_samples + ": stateNumber=" + stateNumber + " lastState=" + lastState + " firstState=" + firstState + " stepSize=" + stepSize);
        lastState = stateNumber;
        return true;
    }

    protected final File file;
    protected final String name;

    private final List<Trace> traces = new ArrayList<Trace>();

    private long burnIn = -1;
    private long firstState = -1;
    private long secondState = -1;
    private long lastState = -1;
    private long stepSize = -1;

    public static class TrimLineReader extends BufferedReader {

        public TrimLineReader(Reader reader) {
            super(reader);
        }

        public String readLine() throws java.io.IOException {
            lineNumber += 1;
            String line = super.readLine();
            if (line != null) return line.trim();
            return null;
        }

        public StringTokenizer tokenizeLine() throws java.io.IOException {
            String line = readLine();
            return getStringTokenizer(line);
        }

        public StringTokenizer getStringTokenizer(String line) {
            if (line == null) return null;
            return new StringTokenizer(line, "\t");
        }

        public int getLineNumber() {
            return lineNumber;
        }

        private int lineNumber = 0;
    }

}

