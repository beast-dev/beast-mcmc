/*
 * ARGTrace.java
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

package dr.evomodel.arg;

import dr.evolution.io.Importer;
import dr.evolution.util.TaxonList;
import dr.util.Identifiable;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.Reader;
import java.util.ArrayList;
import java.util.StringTokenizer;

/**
 * Created by IntelliJ IDEA.
 * User: msuchard
 * Date: Jul 19, 2007
 * Time: 5:30:44 PM
 * To change this template use File | Settings | File Templates.
 */

public class ARGTrace implements Identifiable {

	public ARGTrace() {
	}

	public int getTreeCount(int burnin) {
		int startIndex = (burnin - minState) / stepSize;
		if (startIndex < 0) {
			startIndex = 0;
		}
		return args.size() - startIndex;
	}

	public ARGModel getARG(int index, int burnin) {
		int startIndex = (burnin - minState) / stepSize;
		if (startIndex < 0) {
			startIndex = 0;
		}
		return (ARGModel) args.get(index + startIndex);
	}

	public void add(ARGModel arg) {
		args.add(arg);
	}

	public void setMinimumState(int minState) {
		this.minState = minState;
	}

	public int getMinimumState() {
		return minState;
	}

	public void setStepSize(int stepSize) {
		this.stepSize = stepSize;
	}

	public int getStepSize() {
		return stepSize;
	}

	public int getMaximumState() {
		return (args.size() - 1) * stepSize + minState;
	}

	public String getId() {
		return id;
	}

	public void setId(String id) {
		this.id = id;
	}

	private ArrayList<ARGModel> args = new ArrayList<ARGModel>();

	private int minState;
	private int stepSize;
	private String id;

	/**
	 * Loads the trace for with trees from a reader
	 *
	 * @return the TreeTrace
	 */
	public static ARGTrace loadARGTrace(Reader r) throws IOException, Importer.ImportException {

		BufferedReader reader = new BufferedReader(r);

		ARGTrace trace = new ARGTrace();

		dr.evolution.util.TaxonList taxonList = null;

		int minState = -1;
		int stepSize = 0;

		String line;
		ARGModel nullARG = new ARGModel(null, null, 0, 0);
		ArrayList<String> nameList = new ArrayList<String>();

		taxonList = (TaxonList) nullARG;

		while ((line = reader.readLine()) != null) {

//		String line = reader.readLine();

			line.trim();

			if (line.toUpperCase().startsWith("ARG")) {
				StringTokenizer st = new StringTokenizer(line, "=");
				nameList.add(st.nextToken());
				ARGModel arg = nullARG.fromGraphStringCompressed(st.nextToken().trim());
				trace.add(arg);

			}

		}

		if (nameList.size() < 2) {
			throw new Importer.ImportException("Less than two ARGs in the trace file");
		}

		minState = getStateNumber(nameList.get(0));
		stepSize = getStateNumber(nameList.get(1)) - minState;

//		if (line.toUpperCase().startsWith("#NEXUS")) {
//			NexusImporter importer = new NexusImporter(reader);
//			Tree [] trees = importer.importTrees(null);
//
//			if (trees.length < 2) {
//				throw new Importer.ImportException("Less than two trees in the trace file");
//			}
//
//			String id1 = trees[0].getId();
//            String id2 = trees[1].getId();
//
//            minState = getStateNumber(id1);
//            stepSize = getStateNumber(id2) - minState;
//
//			for (int i = 0; i < trees.length; i++) {
//				args.add(args[i]);
//			}
//		} else {
//			NewickImporter importer = new NewickImporter(reader);
//
//			while (true) {
//
//				int state = 0;
//				Tree tree;
//
//				try {
//					state = importer.readInteger();
//					tree = importer.importTree(taxonList);
//
//					if (taxonList == null) {
//						// The first tree becomes the taxon list. This means
//						// that all subsequent trees will look up their taxa
//						// in that taxon list rather than creating their own
//						// duplicitous ones.
//						taxonList = tree;
//					}
//				} catch (Importer.ImportException ie) {
//					System.out.println("Error reading tree for state " + state);
//					throw ie;
//				} catch (EOFException e) {
//					break;
//				}
//
//				if (minState == -1) {
//					minState = state;
//				} else if (stepSize == 0) {
//					stepSize = state - minState;
//				}
//
//				trace.add(tree);
//			}
//		}

		trace.setMinimumState(minState);
		trace.setStepSize(stepSize);

		return trace;
	}

	private final static int getStateNumber(String id) throws Importer.ImportException {
		try {
			return Integer.parseInt(id.substring(id.indexOf('_') + 1).trim());
		} catch (NumberFormatException nfe) {
			throw new Importer.ImportException("Bad state number in tree label '" + id + "', the state must be preceded by an underscore(_).");
		}
	}
}
