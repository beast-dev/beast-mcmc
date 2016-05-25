/*
 * TracesReader.java
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

import java.io.Reader;
import java.util.ArrayList;
import java.util.List;
import java.util.StringTokenizer;

/**
 * @author Alexei Drummond
 */
public class TracesReader {

    public TracesReader() {
    }

    public void addTracesListener(TracesListener listener) {
        listeners.add(listener);
    }

    public void removeTracesListener(TracesListener listener) {
        listeners.remove(listener);
    }

    /**
     * Reads all the traces in a file, but does not store in memory.
     *
     * @param r the reader to read traces from
     * @throws TraceException      when trace contents is not valid
     * @throws java.io.IOException low level problems with file
     */
    public void readTraces(Reader r) throws TraceException, java.io.IOException {

        LogFileTraces.TrimLineReader reader = new LogFileTraces.TrimLineReader(r);

        // Read through to first token
        StringTokenizer tokens = reader.tokenizeLine();

        // read over empty lines
        while (!tokens.hasMoreTokens()) {
            tokens = reader.tokenizeLine();
        }

        // skip the first column which should be the state number
        String token = tokens.nextToken();

        // lines starting with [ are ignored, assuming comments in MrBayes file
        // lines starting with # are ignored, assuming comments in Migrate or BEAST file
        while (token.startsWith("[") || token.startsWith("#")) {

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
        }
        for (TracesListener listener : listeners) {
            listener.traceNames(labels);
        }

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

            } catch (NumberFormatException nfe) {
                throw new TraceException("State " + state + ":Expected real value in column " + reader.getLineNumber());
            }

            double[] values = new double[labels.length];
            for (int i = 0; i < values.length; i++) {
                if (tokens.hasMoreTokens()) {

                    try {
                        values[i] = Double.parseDouble(tokens.nextToken());
                    } catch (NumberFormatException nfe) {
                        throw new TraceException("State " + state + ": Expected real value in column " + (i + 1) +
                                " (Line " + reader.getLineNumber() + ")");
                    }
                } else {
                    throw new TraceException("State " + state + ": missing values at line " + reader.getLineNumber());
                }
            }

            for (TracesListener listener : listeners) {
                listener.traceRow(state, values);
            }
            tokens = reader.tokenizeLine();
        }
    }

    private List<TracesListener> listeners = new ArrayList<TracesListener>();
}
