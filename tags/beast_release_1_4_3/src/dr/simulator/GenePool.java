/*
 * GenePool.java
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

package dr.simulator;

import dr.simulator.fitnessfunctions.FitnessFunction;
import dr.simulator.mutators.Mutator;
import dr.simulator.datatypes.*;

import java.util.*;

/**
 * @author rambaut
 *         Date: Apr 22, 2005
 *         Time: 9:12:27 AM
 */
public class GenePool {

    public GenePool(int genomeLength, DataType dataType, FitnessFunction fitnessFunction) {
        this.genomeLength = genomeLength;
        this.stateSize = dataType.getStateCount();
	    this.dataType = dataType;
        this.fitnessFunction = fitnessFunction;

        stateFrequencies = new int[genomeLength][];
        for (int i = 0; i < genomeLength; i++) {
            stateFrequencies[i] = new int[stateSize];
        }

        Mutation.initialize(genomeLength, stateSize);
    }

	public Genome createGenome(byte[] sequence) {
        Genome newGenome = new Genome(sequence);

        newGenome.setLogFitness(fitnessFunction.getLogFitness(sequence));

        uniqueGenomeCount++;

        return newGenome;
	}

    /**
     * Replicates a genome with mutations introduced by the given mutator. If no mutations are
     * required then the original genome is returned but with the frequency incremented.
     *
     * @param mutator the mutator object
     * @return the new replicated genome
     */
	public Genome replicateGenome(Genome genome, Mutator mutator) {

        Mutation[] mutations = mutator.mutate(genome);
        if (mutations.length > 0) {
            Genome newGenome;

            if (unusedGenomes.size() > 0) {
                newGenome = (Genome)unusedGenomes.removeFirst();
            } else {
		        newGenome = new Genome();
	            genomes.add(newGenome);
            }

	        newGenome.duplicate(genome);
            newGenome.setFrequency(1);
            newGenome.applyMutations(mutations, fitnessFunction);

            uniqueGenomeCount++;

            return newGenome;
        } else {
            genome.setFrequency(genome.getFrequency() + 1);

            return genome;
        }
	}

    public void killGenome(Genome genome) {

        int frequency = genome.getFrequency() - 1;
	    genome.setFrequency(frequency);
        if (frequency == 0) {
            genome.setFrequency(0);
	        genome.clearMutations();
            unusedGenomes.add(genome);
            uniqueGenomeCount--;
        }
    }

    public void resetStateFrequencies(int populationSize) {
        byte[] masterSequence = Genome.getMasterSequence();

        for (int i = 0; i < masterSequence.length; i++) {
            stateFrequencies[i] = new int[stateSize];
            for (int j = 0; j < stateSize; j++) {
                stateFrequencies[i][j] = 0;
            }
            stateFrequencies[i][masterSequence[i]] = populationSize;
        }
    }

    public void updateStateFrequencies(Genome genome) {
        byte[] masterSequence = Genome.getMasterSequence();

        Mutation[] mutations = genome.getMutations();
		for (int i = 0; i < mutations.length; i++) {
            stateFrequencies[mutations[i].position][mutations[i].state] += 1;
            stateFrequencies[mutations[i].position][masterSequence[mutations[i].position]] -= 1;
		}
    }

    public byte[] getConsensusSequence() {
        byte[] sequence = new byte[genomeLength];

        for (int i = 0; i < genomeLength; i++) {
            sequence[i] = 0;
            int freq = stateFrequencies[i][0];
            for (byte j = 1; j < stateSize; j++) {
                if (freq < stateFrequencies[i][j]) {
                    freq = stateFrequencies[i][j];
                    sequence[i] = j;
                }
            }
        }

        return sequence;
    }

    public void changeMasterSequence(byte[] masterSequence) {
        Iterator iter = genomes.iterator();
        while (iter.hasNext()) {
            Genome genome = (Genome)iter.next();
            genome.changeMasterSequence(masterSequence);
        }
        Genome.setMasterSequence(masterSequence, fitnessFunction.getLogFitness(masterSequence));
    }

	/**
	 * Returns a string containing the DNA sequence of the genome
	 *
	 * @return the DNA sequence
	 */
	public String getDNASequenceString(Genome genome) {
		StringBuffer buffer = new StringBuffer();

		if (dataType instanceof Nucleotides) {
			for (int i = 0; i < genomeLength; i++) {
				buffer.append(dataType.getChar(genome.getState(i)));
			}
		} else if (dataType instanceof Codons) {
			for (int i = 0; i < genomeLength; i++) {
				buffer.append(dataType.getTriplet(genome.getState(i)));
			}
		} else {
			throw new IllegalArgumentException("The genome's data type is not able to be converted to nucleotides");
		}

		return buffer.toString();
	}

	/**
	 * Returns a string containing a translated amino acid sequence of the genome
	 *
	 * @return the amino acid sequence
	 */
	public String getAminoAcidSequenceString(Genome genome) {
		StringBuffer buffer = new StringBuffer();

		if (dataType instanceof AminoAcids) {
			for (int i = 0; i < genomeLength; i++) {
				buffer.append(dataType.getChar(genome.getState(i)));
			}
		} else if (dataType instanceof Codons) {
			for (int i = 0; i < genomeLength; i++) {
				buffer.append(((Codons)dataType).getGeneticCode().getAminoAcidChar(genome.getState(i)));
			}
		} else {
			throw new IllegalArgumentException("The genome's data type is not able to be converted to amino acids");
		}

		return buffer.toString();
	}

    /**
     * Returns the hamming distance (absolute number of differences) between two virusexplorer
     * @param genome1 a genome
     * @param genome2 another genome
     * @return the hamming distance
     */
    public int hammingDistance(Genome genome1, Genome genome2) {

        int distance = 0;
        for (int i = 0; i < genome1.getLength(); i++) {
            if (genome1.getState(i) != genome2.getState(i)) {
                distance += 1;
            }
        }
        return distance;
    }

    public FitnessFunction getFitnessFunction() {
        return fitnessFunction;
    }

    public int getGenomeLength() {
        return genomeLength;
    }

	public List getGenomes() {
	    return genomes;
	}
	
    public int getUniqueGenomeCount() {
        return uniqueGenomeCount;
    }

    public int getUnusedGenomeCount() {
        return unusedGenomes.size();
    }

    public int getStateSize() {
        return stateSize;
    }

    private final FitnessFunction fitnessFunction;

    private final int genomeLength;
	private final DataType dataType;
    private final int stateSize;

    private int uniqueGenomeCount = 0;

    private final int[][] stateFrequencies;

    private final LinkedList genomes = new LinkedList();
    private final LinkedList unusedGenomes = new LinkedList();

}
