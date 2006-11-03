/*
 * Traces.java
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

package dr.app.coalgen;

import java.io.BufferedReader;
import java.io.Reader;
import java.util.ArrayList;
import java.util.StringTokenizer;

/**
 * A class that stores a set of traces from a single chain
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: Traces.java,v 1.2 2004/10/01 22:40:02 alexei Exp $
 */
 
public class Traces implements TraceList {

	public Traces(String name) {
		this.name = name;
	}
	
	/** @return the number of traces in this traceset */
	public int getTraceCount() { return traces.size(); }
	
	/** @return the index of the trace with the given name */
	public int getTraceIndex(String name) {
		for (int i = 0; i < traces.size(); i++) {
			Trace trace = getTrace(i);
			if (name.equals(trace.getName())) {
				return i;
			}
		}
		return -1;
	}
	     
	/** @return the name of the trace with the given index */
	public String getTraceName(int index) {
		return getTrace(index).getName();
	}
	
	/** @return the name of this traceset */
	public String getName() { return name; }

	/** @return the last state in the chain */
	public int getMaxState() {
		return lastState;
	}
	
	/** @return the number of states excluding the burnin */
	public int getStateCount() {
		return ((lastState - getBurnIn()) / stepSize) + 1;
	}
	
	/** @return the size of the step between states */
	public int getStepSize() {
		return stepSize;
	}
	
	public int getBurnIn() { 
		if (burnIn < 0) {
			int n = lastState / stepSize;
			
			setBurnIn((n / 10) * stepSize);
		}
		return burnIn;
	}
	
	public void setBurnIn(int burnIn) { 
		this.burnIn = burnIn;
	}

	public double getStateValue(int trace, int index) {
		return getTrace(trace).getValue(index + (burnIn / stepSize));
	}
		
	public void getValues(int index, double[] destination) {
		 getTrace(index).getValues((burnIn / stepSize), destination, 0);
	}
		
	public void getValues(int index, double[] destination, int offset) {
		 getTrace(index).getValues((burnIn / stepSize), destination, offset);
	}
		
	/**
	 * Loads all the traces in a file
	 */
	public static Traces loadTraces(Reader r, String name) throws TraceException, java.io.IOException {
		
		BufferedReader reader;
		if (r instanceof BufferedReader) {
			reader = (BufferedReader)r;
		} else {
			reader = new BufferedReader(r);
		}
		
		Traces traces = new Traces(name);
		
		// Read through to first token
		String line = reader.readLine().trim();
		StringTokenizer tokens = new StringTokenizer(line,"\t");
		while (!tokens.hasMoreTokens()) {
			line = reader.readLine().trim();
			tokens = new StringTokenizer(line,"\t");
		}
		
		// skip the first column which should be the state number
		String token = tokens.nextToken();
		while (token.startsWith("[")) { // Probably a MrBayes trace file starting with a comment
			line = reader.readLine().trim();
		 	tokens = new StringTokenizer(line,"\t");
			while (!tokens.hasMoreTokens()) {
				 line = reader.readLine().trim();
				 tokens = new StringTokenizer(line,"\t");
			}
			token = tokens.nextToken();
		}
		
		while (tokens.hasMoreTokens()) {
			String label = tokens.nextToken();
			traces.addTrace(label);
		}
		
		int statCount = traces.getTraceCount();

		line = reader.readLine().trim();
	
		while (line != null && !line.equals("")) {
		
			tokens = new StringTokenizer(line);
			if (!tokens.hasMoreTokens()) break;
			
			String stateString = tokens.nextToken();
			int state = 0;
			try {
				state = Integer.parseInt(stateString);
				
				// MrBayes puts 1 as the first state, BEAST puts 0
				// In order to get the same gap between subsequent samples,
				// we force this to 0.
				if (state == 1) state = 0;
				traces.addState(state);
	
			} catch (NumberFormatException nfe) {
				break;
			}
				
			for (int i = 0; i < statCount; i++) {
				if (tokens.hasMoreTokens()) {
				
					try {
						traces.addValue(i, Double.parseDouble(tokens.nextToken()));
					} catch (NumberFormatException nfe) {
						throw new TraceException("State " + state + ": Expected real value in column " + String.valueOf(i + 2));
					}
				} else {
					throw new TraceException("State " + state + ": Expected real value in column " + String.valueOf(i + 2));
				}
				
			}
			
			line = reader.readLine();
			if (line != null) line = line.trim();
		}
		 
		return traces;
	}

	
	//************************************************************************
	// private methods
	//************************************************************************
	
	// These methods are used by the load function, above
	
	/**
	 * Add a trace for a statistic of the given name
	 */
	private void addTrace(String name) {
		traces.add(new Trace(name));
	}
	
	/**
	 * Add a state number for these traces. This should be
	 * called before adding values for each trace. The spacing 
	 * between stateNumbers should remain constant.
	 */
	private void addState(int stateNumber) throws TraceException {
		if (firstState < 0) {
			firstState = stateNumber;
		} else if (stepSize < 0) {
			stepSize = stateNumber - firstState;
		} else {
			int step = stateNumber - lastState;
			if (step != stepSize) {
				throw new TraceException("State step sizes are not constant");
			}
		}
		lastState = stateNumber;
	}
	
	/**
	 * add a value for the ith trace
	 */
	private void addValue(int i, double value) {
		getTrace(i).add(value);
	}
	
	/**
	 * @return the trace for a given index
	 */
	private Trace getTrace(int index) {
		return (Trace)traces.get(index);
	}
	
	private ArrayList traces = new ArrayList();
	
	private String name;
	
	private int burnIn = -1;
	private int firstState = -1;
	private int lastState = -1;
	private int stepSize = -1;
}