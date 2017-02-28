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
import java.util.TreeMap;

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
        return (getBurnIn() / stepSize);
    }

    /**
     * @return the size of the step between states
     */
    public int getStepSize() {
        return stepSize;
    }

    public int getBurnIn() {
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

    public void setBurnIn(long burnin2) {
        this.burnIn = (int) burnin2;
        for (Trace trace : traces) {
            trace.setTraceStatistics(null);
        }
    }

    public double getStateValue(int trace, int index) {
        return (Double) getTrace(trace).getValue(index + (burnIn / stepSize));
    }

    /**
     * Read several consecutive values of one state into a destination array
     *
     * @param nState      State index number
     * @param destination array to store result
     * @param offset      first trace index
     */
    public void getStateValues(int nState, double[] destination, int offset) {
        final int index1 = nState + (burnIn / stepSize);
        for (int k = 0; k < destination.length; ++k) {
            destination[k] = (Double) getTrace(k + offset).getValue(index1);
        }
    }

    public List getValues(int index, int fromIndex, int toIndex) {
        List newList = null;
        try {
            newList = getTrace(index).getValues(fromIndex, toIndex);
        } catch (Exception e) {
            System.err.println("getValues error: trace index = " + index);
        }
        return newList;
    }

    public List getValues(int index) {
        return this.getValues(index, getBurninStateCount(), getTrace(index).getValueCount());
    }

    public List getBurninValues(int index) {
        return this.getValues(index, 0, getBurninStateCount());
    }

    public void loadTraces() throws TraceException, IOException {
        FileReader reader = new FileReader(file);
        loadTraces(reader);
        reader.close();
    }

    /**
     * Walter: Please comment what the extra arguments mean
     *
     * @param r
     * @throws TraceException
     * @throws java.io.IOException
     */
    public void loadTraces(Reader r) throws TraceException, java.io.IOException {

        TrimLineReader reader = new LogFileTraces.TrimLineReader(r);

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
            readTraceType(token, tokens); // using # to define type
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
            addTraceAndType(labels[i]);
        }


        int traceCount = getTraceCount();

        int num_samples = 0;

        tokens = reader.tokenizeLine();
        while (tokens != null && tokens.hasMoreTokens()) {

            String stateString = tokens.nextToken();
            long state = 0;

            try {
                try {
                    // Changed this to parseDouble because LAMARC uses scientific notation for the state number
                    state = (long) Double.parseDouble(stateString);
                } catch (NumberFormatException nfe) {
                    throw new TraceException("Unable to parse state number in column 1 (Line " + reader.getLineNumber() + ")");
                }

                if (num_samples < 1) {
                    // MrBayes puts 1 as the first state, BEAST puts 0
                    // In order to get the same gap between subsequent samples,
                    // we force this to 0.
                    if (state == 1) state = 0;
                }
                num_samples += 1;

                if (!addState(state, num_samples)) {
                    throw new TraceException("State " + state + " is not consistent with previous spacing (Line " + reader.getLineNumber() + ")");
                }

            } catch (NumberFormatException nfe) {
                throw new TraceException("State " + state + ":Expected real value in column " + reader.getLineNumber());
            }

            for (int i = 0; i < traceCount; i++) {
                if (tokens.hasMoreTokens()) {
                    String value = tokens.nextToken();

                    if (state == 0) assignTraceTypeAccordingValue(value);

                    try {
//                        values[i] = Double.parseDouble(tokens.nextToken());
                        addParsedValue(i, value);
                    } catch (NumberFormatException nfe) {
                        throw new TraceException("State " + state + ": Expected correct number type (Double, Integer or String) in column "
                                + (i + 1) + " (Line " + reader.getLineNumber() + ")");
                    }

                } else {
                    throw new TraceException("State " + state + ": missing values at line " + reader.getLineNumber());
                }
            }

            tokens = reader.tokenizeLine();
        }

        burnIn = (int) (0.1 * lastState);
    }

    /**
     * add a value for the n'th trace
     *
     * @param nTrace trace index
     * @param value  next value
     */
    private void addParsedValue(int nTrace, String value) {
        String name = getTraceName(nTrace);
//        System.out.println(thisTrace.getTraceType() + "   " + value);
        if (tracesType.get(name).isNumber()) {
            Double v = Double.parseDouble(value);
            getTrace(nTrace).add(v);

        } else  {
            getTrace(nTrace).add(value);
        }
    }

    private void assignTraceTypeAccordingValue(String value) {
        //todo
    }

    /**
     * todo I would have thought this would read the type across the header for each trace
     * @param firstToken
     * @param tokens
     */
    private void readTraceType(String firstToken, StringTokenizer tokens) {
        if (tokens.hasMoreTokens()) {
            String token; //= tokens.nextToken();
            if (firstToken.toLowerCase().contains(TraceType.ORDINAL.toString())) {
                while (tokens.hasMoreTokens()) {
                    token = tokens.nextToken();
                    tracesType.put(token, TraceType.ORDINAL);
                }
            } else if (firstToken.toLowerCase().contains(TraceType.CATEGORICAL.toString())) {
                while (tokens.hasMoreTokens()) {
                    token = tokens.nextToken();
                    tracesType.put(token, TraceType.CATEGORICAL);
                }
            }
        }
    }

//    public Trace<?> assignTraceType(String name, int numberOfLines) {
//        Trace<?> trace = null;
//        if (tracesType != null) {
//            if (tracesType.get(name) == TraceFactory.TraceType.INTEGER) {
//                trace = new DiscreteTrace(name, numberOfLines);
////                trace.setTraceType(TraceType.INTEGER);
//            } else if (tracesType.get(name) == TraceFactory.TraceType.STRING) {
//                trace = new CategoryTrace(name, numberOfLines);
////                trace.setTraceType(TraceType.STRING);
//            }
//        } else {
//            trace = new ContinuousTrace(name, numberOfLines); // default DOUBLE
//        }
//        return trace;
//    }


    //************************************************************************
    // private methods
    //************************************************************************

    // These methods are used by the load function, above

    /**
     * Add a trace for a statistic of the given name
     *
     * @param name trace name
     */
    private void addTraceAndType(String name) {
        if (tracesType.get(name) == null) {
            traces.add(createTrace(name, TraceType.REAL));
            tracesType.put(name, TraceType.REAL);
        } else {
            traces.add(createTrace(name, tracesType.get(name)));
        }
    }

    private Trace createTrace(String name, TraceType traceType) {
        if (traceType.isNumber()) {
            return new Trace<Double>(name, traceType);
        } else {
            return new Trace<String>(name, TraceType.ORDINAL);
        }
    }

    // TODO get rid of generic to make things easy
    // TODO change to String only, and parse to double, int or string in getValues according to trace type
    public void changeTraceType(int id, TraceType newType) throws TraceException {
        if (id >= getTraceCount() || id < 0) throw new TraceException("trace id is invaild " + id);
        Trace trace = traces.get(id);
        if (trace.getTraceType() != newType) {
            Trace newTrace = null;
            try {
                if (trace.getTraceType() == TraceType.CATEGORICAL) {
                    if (newType == TraceType.REAL) {
                        newTrace = createTrace(trace.getName(), newType);
                    }
                    for (int i = 0; i < trace.getValueCount(); i++) { // String => Double
                        newTrace.add(Double.parseDouble(trace.getValue(i).toString()));
                    }
                } else if (newType == TraceType.CATEGORICAL) {
                    newTrace = createTrace(trace.getName(), TraceType.CATEGORICAL);
                    for (int i = 0; i < trace.getValueCount(); i++) { // Double => String
                        newTrace.add(trace.getValue(i).toString());
                    }
                } else {
                    newTrace = createTrace(trace.getName(), newType); // not need to copy values, because they are both Double
                }
            } catch (Exception e) {
                throw new TraceException("Type change is failed, when parsing " + trace.getTraceType()
                        + " to " + newType + " in trace " + trace.getName());
            }

            if (trace.getTraceType() == TraceType.CATEGORICAL || newType == TraceType.CATEGORICAL) {
                if (newTrace.getValueCount() != trace.getValueCount())
                    throw new TraceException("Type change is failed, because values size is different after copy !");

                traces.set(id, newTrace);
            } else {
                trace.setTraceType(newType);
            }
        }
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
    private boolean addState(long stateNumber, int num_samples) {
        if (firstState < 0) { // it can use num_samples==1 to replace firstState < 0
            firstState = stateNumber;
        }
        // because of forcing MrBayes 1st state=0, initialising from the 1st two rows are not correct
        if (num_samples > 2) {
            if (num_samples == 3)
                stepSize = (int) (stateNumber - lastState);
            // validate if step sizes are same, it needs 2 rows to initiate properly
            int step = (int) (stateNumber - lastState);
            if (step != stepSize) {
                System.out.println("step: " + step + " != " + stepSize);
                return false;
            }
        }
        System.out.println(num_samples + ": stateNumber=" + stateNumber + " lastState=" + lastState + " firstState=" + firstState + " stepSize=" + stepSize);
        lastState = stateNumber;
        return true;
    }

    protected final File file;
    protected final String name;

    private final List<Trace> traces = new ArrayList<Trace>();

    public void addTrace(String newTName, int i) {
        TraceCustomized tc = new TraceCustomized(newTName);
        tc.addValues(traces.get(i)); // only Double
        traces.add(tc);
        tracesType.put(newTName, TraceType.REAL);
    }

    // tracesType only save INTEGER and STRING, and only use during loading files
    private TreeMap<String, TraceType> tracesType = new TreeMap<String, TraceType>();

    private int burnIn = -1;
    private long firstState = -1;
    private long lastState = -1;
    private int stepSize = -1;

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
            if (line == null) return null;
            return new StringTokenizer(line, "\t");
        }

        public int getLineNumber() {
            return lineNumber;
        }

        private int lineNumber = 0;
    }

//    public class D extends LogFileTraces implements TraceList.D {
//
//        public D(String name, File file) {
//            super(name, file);
//        }
//
//        public Double[] getValues(int index, int length) {
//            return this.getValues(index, length, 0);
//        }
//
//        public Double[] getValues(int index, int length, int offset) {
//            Double[] destination = null;
//            try {
//                destination = ((Trace.D) getTrace(index)).getValues(length, getBurninStateCount(), offset, selected);
//            } catch (Exception e) {
//                System.err.println("getValues error: trace index = " + index);
//            }
//            return destination;
//        }
//
//        public Double[] getBurninValues(int index, int length) {
//            Double[] destination = null;
//            try {
//                destination = (Double[]) getTrace(index).getValues(length, 0, 0, getBurninStateCount(), selected);
//            } catch (Exception e) {
//                System.err.println("getValues error: trace index = " + index);
//            }
//            return destination;
//        }
//    }
}

