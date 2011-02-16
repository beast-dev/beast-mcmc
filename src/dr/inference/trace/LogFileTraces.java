/*
 * LogFileTraces.java
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
    public int getMaxState() {
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
        return ((lastState - firstState) / stepSize) - (getBurnIn() / stepSize) + 1;
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

    public void setBurnIn(int burnIn) {
        this.burnIn = burnIn;
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

//    public Object[] createValues(int index, int length) {
//        Trace trace = getTrace(index);
//        return getTrace(index).createValues((burnIn / stepSize), length);
//    }

    public <T> void getValues(int index, T[] destination) {
        try {
            getTrace(index).getValues(getBurninStateCount(), destination, 0);
        } catch (Exception e) {
            System.err.println("getValues error: trace index = " + index);
        }
    }

    public <T> void getValues(int index, T[] destination, int offset) {
        ((Trace<T>) getTrace(index)).getValues(getBurninStateCount(), destination, offset);
    }

    public <T> void getBurninValues(int index, T[] destination) {
        try {
            ((Trace<T>) getTrace(index)).getValues(0, getBurninStateCount(), destination, 0);
        } catch (Exception e) {
            System.err.println("getValues error: trace index = " + index);
        }
    }


    public List getValues(int index, int fromIndex, int toIndex) {
        List newList = null;
        try {
            newList = getTrace(index).getValues(fromIndex, toIndex, selected);
        } catch (Exception e) {
            System.err.println("getValues error: trace index = " + index);
        }
        return newList;
    }

    public List getValues(int index) {
        return this.getValues(index, getBurninStateCount(), getTrace(index).getValuesSize());
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

        boolean firstState = true;

        tokens = reader.tokenizeLine();
        while (tokens != null && tokens.hasMoreTokens()) {

            String stateString = tokens.nextToken();
            int state = 0;

            try {
                try {
                    // Changed this to parseDouble because LAMARC uses scientific notation for the state number
                    state = (int) Double.parseDouble(stateString);
                } catch (NumberFormatException nfe) {
                    throw new TraceException("Unable to parse state number in column 1 (Line " + reader.getLineNumber() + ")");
                }

                if (firstState) {
                    // MrBayes puts 1 as the first state, BEAST puts 0
                    // In order to get the same gap between subsequent samples,
                    // we force this to 0.
                    if (state == 1) state = 0;
                    firstState = false;
                }

                if (!addState(state)) {
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
        if (tracesType.get(name) == TraceFactory.TraceType.CONTINUOUS) {
            Double v = Double.parseDouble(value);
            getTrace(nTrace).add(v);

        } else if (tracesType.get(name) == TraceFactory.TraceType.INTEGER) {
//             Integer v = Integer.parseInt(value);
            int v = (int) Double.parseDouble(value);
            getTrace(nTrace).add(v);

        } else if (tracesType.get(name) == TraceFactory.TraceType.CATEGORY) {
            getTrace(nTrace).add(value);

        } else {
            throw new RuntimeException("Trace type is not recognized: " + tracesType.get(name));
        }
    }

    private void assignTraceTypeAccordingValue(String value) {
        //todo
    }

    private void readTraceType(String firstToken, StringTokenizer tokens) {
        if (tokens.hasMoreTokens()) {
            String token; //= tokens.nextToken();
            if (firstToken.contains(TraceFactory.TraceType.INTEGER.toString())
                    || firstToken.contains(TraceFactory.TraceType.INTEGER.toString().toUpperCase())) {
                while (tokens.hasMoreTokens()) {
                    token = tokens.nextToken();
                    tracesType.put(token, TraceFactory.TraceType.INTEGER);
                }
            } else if (firstToken.contains(TraceFactory.TraceType.CATEGORY.toString())
                    || firstToken.contains(TraceFactory.TraceType.CATEGORY.toString().toUpperCase())) {
                while (tokens.hasMoreTokens()) {
                    token = tokens.nextToken();
                    tracesType.put(token, TraceFactory.TraceType.CATEGORY);
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
//            } else if (tracesType.get(name) == TraceFactory.TraceType.CATEGORY) {
//                trace = new CategoryTrace(name, numberOfLines);
////                trace.setTraceType(TraceType.CATEGORY);
//            }
//        } else {
//            trace = new ContinuousTrace(name, numberOfLines); // default CONTINUOUS
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
        if (tracesType == null || tracesType.get(name) == null || tracesType.get(name) == TraceFactory.TraceType.CONTINUOUS) {
            traces.add(createTrace(TraceFactory.TraceType.CONTINUOUS, name));
            tracesType.put(name, TraceFactory.TraceType.CONTINUOUS);
        } else {
            traces.add(createTrace(tracesType.get(name), name));
        }
    }

//    private void replaceTrace(String name, TraceFactory.TraceType reloadType, int numberOfLines) throws TraceException {
//        if (reloadType == null) throw new TraceException("trace (" + name + ") type is null");
//        if (numberOfLines < 0) throw new TraceException("numberOfLines cannot be 0");
//
//        for (int i = 0; i < traces.size(); i++) {
//            Trace trace = traces.get(i);
//            // change trace type
//            if (trace.getName().equalsIgnoreCase(name)) {
//                if (trace.getTraceType() != reloadType.getType())
//                    traces.set(i, createTrace(reloadType, name));
//                return;
//            }
//        }
//        throw new TraceException("Cannot find trace : " + name);
//    }

    private Trace createTrace(TraceFactory.TraceType traceType, String name) {
        // System.out.println("create trace (" + name + ") with type " + traceType);
        switch (traceType) {
            case CONTINUOUS:
                return new Trace<Double>(name);
            case INTEGER:
                return new Trace<Integer>(name);
            case CATEGORY:
                return new Trace<String>(name);
        }
        throw new IllegalArgumentException("The trace type " + traceType + " is not recognized.");
    }

    /**
     * Add a state number for these traces. This should be
     * called before adding values for each trace. The spacing
     * between stateNumbers should remain constant.
     *
     * @param stateNumber the state
     * @return false if the state number is inconsistent
     */
    private boolean addState(int stateNumber) {
        if (firstState < 0) {
            firstState = stateNumber;
        } else if (stepSize < 0) {
            stepSize = stateNumber - firstState;
        } else {
            int step = stateNumber - lastState;
            if (step != stepSize) {
                return false;
            }
        }
        lastState = stateNumber;
        return true;
    }

    protected final File file;
    protected final String name;

    private final List<Trace> traces = new ArrayList<Trace>();
    // List traces is added before having types
    // tracesType only save INTEGER and CATEGORY, and only use during loading files
    private TreeMap<String, TraceFactory.TraceType> tracesType = new TreeMap<String, TraceFactory.TraceType>();

    private int burnIn = -1;
    private int firstState = -1;
    private int lastState = -1;
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

