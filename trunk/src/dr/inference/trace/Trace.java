/*
 * Trace.java
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

import dr.util.Identifiable;

import java.io.BufferedReader;
import java.io.Reader;
import java.util.ArrayList;
import java.util.StringTokenizer;

/**
 *
 * @version
 * @author Alexei Drummond
 */
public interface Trace extends Identifiable {
	
	double[] getValues(int burnIn);
	
	/**
	 * Copies the trace values (after discarding burnin) into a destination array at a given start index
	 */
	void copyValues(int burnIn, double[] destination, int startIndex);
	
	int getMinimumState();
	
	int getStepSize();
	
	int getValueCount(int burnin);
	
	class Utils {
		
		public static double[] getStates(Trace trace, int burnin) {
			double[] states = new double[trace.getValueCount(burnin)];
			int state = trace.getMinimumState() + burnin;
			for (int i =0; i < states.length; i++) {
				states[i] = state;
				state += trace.getStepSize();
			}
			return states;
		}
		
		public static int getMaximumState(Trace trace) {
			return (trace.getValueCount(0)-1) * trace.getStepSize() + trace.getMinimumState();
		}
		
		/**
		 * Loads the trace for a single statistic given by name
		 * @return null if statistic not found in trace file
		 */
		public static Trace loadTrace(Reader r, String name) throws java.io.IOException {
			
			BufferedReader reader;
			if (r instanceof BufferedReader) {
				reader = (BufferedReader)r;
			} else {
				reader = new BufferedReader(r);
			}
			
			int statNum = -1;
			
			// Read through to first token
			String line = reader.readLine().trim();
			StringTokenizer tokens = new StringTokenizer(line,"\t");
			while (!tokens.hasMoreTokens()) {
				 line = reader.readLine().trim();
				 tokens = new StringTokenizer(line,"\t");
			}
			
			// skip state token
			tokens.nextToken();
			int statCount = 0;
			while (tokens.hasMoreTokens()) {
				
				if (name.equals(tokens.nextToken())) {
					statNum = statCount;
				}
				statCount++;
			}
			
			if (statNum == -1) return null;
			
			int minState = -1;
			int stepSize = 0;
			ArrayList values = new ArrayList();
			
			line = reader.readLine().trim();
		
			while (line != null && !line.equals("")) {
			
				tokens = new StringTokenizer(line);
				if (!tokens.hasMoreTokens()) break;
				String stateString = tokens.nextToken();
				try {
					int state = Integer.parseInt(stateString);
					if (minState == -1) { 
						minState = state; 
					} else if (stepSize == 0) {
						stepSize = state - minState;
					}
					for (int i = 0; i < statCount; i++) {
						if (i == statNum) {
							values.add(new Double(tokens.nextToken()));
						} else {
							tokens.nextToken();
						}
					}
					line = reader.readLine();
				} catch (NumberFormatException nfe) {
					break;
				}
				if (line != null) line = line.trim();
			}
			Trace trace = new DefaultTrace(name, values, minState, stepSize);
			
			return trace;
		}
	
	
		/**
		 * Loads all the traces in a file
		 */
		public static Trace[] loadTraces(Reader r) throws java.io.IOException {
			
			BufferedReader reader;
			if (r instanceof BufferedReader) {
				reader = (BufferedReader)r;
			} else {
				reader = new BufferedReader(r);
			}
			
			ArrayList names = new ArrayList();
			
			// Read through to first token
			String line = reader.readLine().trim();
			StringTokenizer tokens = new StringTokenizer(line,"\t");
			while (!tokens.hasMoreTokens()) {
				 line = reader.readLine().trim();
				 tokens = new StringTokenizer(line,"\t");
			}
			
			// skip tokens up to state token
			while (tokens.hasMoreTokens() && !tokens.nextToken().toLowerCase().endsWith("state")) {}
			
			while (tokens.hasMoreTokens()) {
				String name = tokens.nextToken();
				if (name.length() != 0) {
					names.add(name);
				}
			}
			
			int statCount = names.size();
			int minState = -1;
			int stepSize = 0;
			ArrayList[] values = new ArrayList[statCount];
			for (int i =0; i < values.length; i++) {
				values[i] = new ArrayList();
			}
			
			line = reader.readLine().trim();
		
			while (line != null && !line.equals("")) {
			
				tokens = new StringTokenizer(line);
				if (!tokens.hasMoreTokens()) break;
				String stateString = tokens.nextToken();
				try {
					int state = Integer.parseInt(stateString);
					if (minState == -1) { 
						minState = state; 
					} else if (stepSize == 0) {
						stepSize = state - minState;
					}
					for (int i = 0; i < statCount; i++) {
						if (tokens.hasMoreTokens()) {
							values[i].add(new Double(tokens.nextToken()));
						} else {
							StringBuffer buf = new StringBuffer("State " + state + ", Expected double after " );
							for (int j = 0; j < Math.min(i, 10); j++) {
								Double d = (Double)values[j].get(values[j].size()-1);
								buf.append(names.get(j)+":"+d.doubleValue()+" ");
							}
							throw new RuntimeException(buf.toString());
						}
					}
					line = reader.readLine();
				} catch (NumberFormatException nfe) {
					break;
				}
				if (line != null) line = line.trim();
			}
			Trace[] traces = new Trace[statCount];
			for (int i = 0; i < traces.length; i++) {
				traces[i] = new DefaultTrace((String)names.get(i), values[i], minState, stepSize);
			}
			return traces;
		}
	}
	
	class DefaultTrace implements Trace {
		
		
		public DefaultTrace(String id, double[] values, int minState, int stepSize) {
			this.values = values;
			this.id = id;
			this.minState = minState;
			this.stepSize = stepSize;
		}
		
		public DefaultTrace(String id, ArrayList valueList, int minState, int stepSize) {
		
			values = new double[valueList.size()];
			for (int i =0; i < values.length; i++) {
				values[i] = ((Double)valueList.get(i)).doubleValue();
			}
			this.id = id;
			this.minState = minState;
			this.stepSize = stepSize;	
		}
	
		public double[] getValues(int burnIn) {
			int startIndex = (burnIn - getMinimumState()) / getStepSize();
			double[] newValues = new double[values.length-startIndex];
			System.arraycopy(values, startIndex, newValues, 0, newValues.length);
			return newValues;
		}
		
		/**
		 * Copies the trace values (after discarding burnin) into a destination array at a given start index
		 */
		public void copyValues(int burnIn, double[] destination, int destIndex) {
			int startIndex = (burnIn - getMinimumState()) / getStepSize();
			if (destination.length - destIndex < (values.length - startIndex)) { 
				throw new RuntimeException("Destination array not big enough"); 
			}
			System.arraycopy(values, startIndex, destination, destIndex, values.length - startIndex);
		}
		
		public int getMinimumState() { return minState; }
		public int getStepSize() { return stepSize; }
		public int getValueCount() { return values.length; }
		
		public int getValueCount(int burnin) { 
			int startIndex = (burnin - getMinimumState()) / getStepSize();
			return values.length - startIndex;
		}

		public String getId() { return id; }
		public void setId(String id) { this.id = id; }
		
		private double[] values;
		private int minState;
		private int stepSize;
		private String id;
	}
}
