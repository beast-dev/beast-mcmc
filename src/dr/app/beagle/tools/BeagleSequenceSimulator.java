/*
 * BeagleSequenceSimulator.java
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

package dr.app.beagle.tools;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import dr.app.bss.Utils;
import dr.evolution.alignment.SimpleAlignment;
import dr.evolution.datatype.DataType;
import dr.evolution.sequence.Sequence;
import dr.evolution.tree.NodeRef;
import dr.evolution.util.Taxon;

/**
 * @author Filip Bielejec
 * @version $Id$
 * 
 */
public class BeagleSequenceSimulator {

	// Constructor fields
	private ArrayList<Partition> partitions;
	private int siteCount;

	// Alignment fields
	public static final int gapFlag = Integer.MAX_VALUE;
	private SimpleAlignment alignment;
	private DataType dataType;
	private boolean fieldsSet = false;
	LinkedHashMap<Integer, LinkedHashMap<NodeRef, int[]>> partitionSequencesMap;
//	private boolean outputAncestralSequences = true;
	
	public BeagleSequenceSimulator(ArrayList<Partition> partitions) {

		this.partitions = partitions;
		this.alignment = new SimpleAlignment();
		alignment.setReportCountStatistics(false);

		partitionSequencesMap = new LinkedHashMap<Integer, LinkedHashMap<NodeRef,int[]>>();
		
		int siteCount = 0;
		int to = 0;
		for (Partition partition : partitions) {

			to = partition.to;
			if (to > siteCount) {
				siteCount = to;
			}

			if (!fieldsSet) {
				
				dataType = partition.getDataType();
//				outputAncestralSequences = partition.isOutputAncestralSequences();
				fieldsSet = true;
//				System.err.println(dataType);
				
			} else {

				if (dataType.getType() != partition.getDataType().getType()) {
					throw new RuntimeException(
							"Partitions must have the same data type.");
				}

			}

		}// END: partitions loop

		this.siteCount = siteCount + 1;
	}// END: Constructor

	public SimpleAlignment simulate(boolean parallel, boolean outputAncestralSequences) {

		try {

			// Executor for threads
			int NTHREDS = 1;
			if (parallel) {
				NTHREDS = Math.min(partitions.size(), Runtime.getRuntime()
						.availableProcessors());
			}

			ExecutorService executor = Executors.newFixedThreadPool(NTHREDS);
			List<Callable<Void>> simulatePartitionCallers = new ArrayList<Callable<Void>>();

			int partitionCount = 0;
			for (Partition partition : partitions) {

				partition.setPartitionNumber(partitionCount);
                partition.setOutputAncestralSequences(outputAncestralSequences);
				
				simulatePartitionCallers.add(new SimulatePartitionCallable(
						partition
//						, partitionCount
						));
				partitionCount++;

			}// END: partitions loop

			executor.invokeAll(simulatePartitionCallers);

			// Wait until all threads are finished
			executor.shutdown();
			while (!executor.isTerminated()) {
			}

			alignment = compileAlignment();

		} catch (Exception e) {
			e.printStackTrace();
		}// END: try-catch block

		return alignment;
	}// END: simulate

	private class SimulatePartitionCallable implements Callable<Void> {

		private Partition partition;
//        private int partitionNumber;
		
		private SimulatePartitionCallable(Partition partition
//				, int partitionNumber
				) {
			this.partition = partition;
//			this.partitionNumber = partitionNumber;
		}// END: Constructor

		public Void call() {

			try {

				partition.simulatePartition();
//                partitionSequencesMap.put(partitionNumber, partition.getSequenceMap());
                
			} catch (Exception e) {
				Utils.handleException(e);
			}

			return null;
		}// END: call

	}// END: SimulatePartitionCallable class

	private SimpleAlignment compileAlignment() {

		SimpleAlignment simpleAlignment = new SimpleAlignment();
		simpleAlignment.setReportCountStatistics(false);
		simpleAlignment.setDataType(dataType);
		
		LinkedHashMap<Taxon, int[]> alignmentMap = new LinkedHashMap<Taxon, int[]>();

		// compile the alignment
		for (Partition partition : partitions) {

			Map<Taxon, int[]> sequenceMap = partition.getTaxonSequencesMap();
			Iterator<Entry<Taxon, int[]>> iterator = sequenceMap.entrySet()
					.iterator();

			while (iterator.hasNext()) {

				Entry<Taxon, int[]> pairs = (Entry<Taxon, int[]>) iterator
						.next();

				Taxon taxon = pairs.getKey();
				int[] partitionSequence = pairs.getValue();

				if (alignmentMap.containsKey(taxon)) {

					int j = 0;
					for (int i = partition.from; i <= partition.to; i += partition.every) {

						alignmentMap.get(taxon)[i] = partitionSequence[j];
						j++;

					}// END: i loop

				} else {

					int[] sequence = new int[siteCount];
					// dirty solution for gaps when taxa between the tree
					// topologies don't match
					Arrays.fill(sequence, gapFlag);

					int j = 0;
					for (int i = partition.from; i <= partition.to; i += partition.every) {

						sequence[i] = partitionSequence[j];
						j++;

					}// END: i loop

					alignmentMap.put(taxon, sequence);
				}// END: key check

			}// END: iterate seqMap

		}// END: partitions loop

		Iterator<Entry<Taxon, int[]>> iterator = alignmentMap.entrySet()
				.iterator();
		while (iterator.hasNext()) {

			Entry<Taxon, int[]> pairs = (Entry<Taxon, int[]>) iterator.next();
			Taxon taxon = (Taxon) pairs.getKey();
			int[] intSequence = (int[]) pairs.getValue();
			
			Sequence sequence = Utils.intArray2Sequence(taxon, //
					intSequence, //
					gapFlag, //
					dataType
			);
			
//			sequence.setDataType(dataType);
			
			simpleAlignment.addSequence(sequence);

			iterator.remove();

		}// END: while has next

		return simpleAlignment;
	}// END: compileAlignment

	public LinkedHashMap<Integer, LinkedHashMap<NodeRef, int[]>> getPartitionSequencesMap() {
		return partitionSequencesMap;
	}//END: getPartitionSequencesMap
	
//	private Sequence intArray2Sequence(Taxon taxon, int[] seq, int gapFlag) {
//
//		StringBuilder sSeq = new StringBuilder();
//
//		if (dataType instanceof Codons) {
//
//			for (int i = 0; i < siteCount; i++) {
//
//				int state = seq[i];
//
//				if (state == gapFlag) {
//					sSeq.append(dataType.getTriplet(dataType.getGapState()));
//				} else {
//					sSeq.append(dataType.getTriplet(seq[i]));
//				}// END: gap check
//
//			}// END: replications loop
//
//		} else {
//
//			for (int i = 0; i < siteCount; i++) {
//
//				int state = seq[i];
//
//				if (state == gapFlag) {
//					sSeq.append(dataType.getCode(dataType.getGapState()));
//				} else {
//					sSeq.append(dataType.getCode(seq[i]));
//				}// END: gap check
//
//			}// END: replications loop
//
//		}// END: dataType check
//
//		return new Sequence(taxon, sSeq.toString());
//	}// END: intArray2Sequence

} // END: class
