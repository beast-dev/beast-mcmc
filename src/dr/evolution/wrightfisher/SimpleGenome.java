/*
 * SimpleGenome.java
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

package dr.evolution.wrightfisher;

import dr.evolution.datatype.GeneticCode;

public class SimpleGenome implements Genome {

	public SimpleGenome(int size, FitnessFunction fitnessFunction, boolean initializeToFittest) {
		this.sequence = new byte[size];
		if (initializeToFittest) fitnessFunction.initializeToFittest(sequence);
		this.fitness = 1.0;
		this.parent = null;
	}

	private SimpleGenome(byte[] sequence, double fitness, Genome parent, Mutation[] mutations) {
		this.sequence = sequence;
		this.fitness = fitness;
		this.parent = parent;
		this.mutations = mutations;
	}

	public Genome replicate(Mutator mutator, FitnessFunction fitnessFunction) {

		byte[] childSequence = new byte[sequence.length];

		Mutation[] mutations = mutator.mutate(sequence, childSequence);

		double newFitness = fitnessFunction.getFitness(childSequence);

		return new SimpleGenome(childSequence, newFitness, this, mutations);
	}

	public byte[] getSequence() {
		return sequence;
	}

	public int hammingDistance(Genome other) {

		int distance = 0;
		byte[] otherSequence = other.getSequence();
		for (int i = 0; i < sequence.length; i++) {
			if (sequence[i] != otherSequence[i]) distance += 1;
		}
		return distance;
	}


	public String getDNASequenceString() {
		StringBuffer buffer = new StringBuffer();
		for (int i = 0; i < sequence.length; i++) {
			switch (sequence[i]) {
				case 0: buffer.append('A'); break;
				case 1: buffer.append('C'); break;
				case 2: buffer.append('G'); break;
				case 3: buffer.append('T'); break;
				default: buffer.append('N'); break;
			}
		}
		return buffer.toString();
	}

	public String getAminoAcidSequenceString() {
		StringBuffer buffer = new StringBuffer();
		for (int i = 0; i < sequence.length; i+=3) {
			int codonState = (sequence[i] * 16) + (sequence[i+1] * 4) + sequence[i+2];
			//int codonState = codons.getState(sequence[i], sequence[i+1], sequence[i+2]);

			char aminoAcid = GeneticCode.UNIVERSAL.getAminoAcidChar(codonState);

			buffer.append(aminoAcid);
		}
		return buffer.toString();
	}

	public final void mark() {
		marks += 1;
		if (parent != null) parent.mark();
	}

	public final void unmark() {
		marks = 0;
		if (parent != null) parent.unmark();
	}

	public final int getMarks() { return marks; }

    public final Genome getParent() { return parent; }

	public final double getFitness() { return fitness; }
    public final Mutation[] getMutations() { return mutations; }

	public final int getGenomeLength() { return sequence.length; }

	byte[] sequence;
	double fitness;
	Mutation[] mutations = null;
	Genome parent;
	int marks = 0;
}