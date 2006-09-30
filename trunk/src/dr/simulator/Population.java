/*
 * Population.java
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

import dr.simulator.fitnessfunctions.*;
import dr.simulator.mutators.Mutator;
import dr.simulator.mutators.CodonMutator;
import dr.simulator.selectors.Selector;
import dr.simulator.selectors.NaturalSelector;
import dr.simulator.datatypes.GeneticCode;
import dr.simulator.datatypes.Codons;


/**
 * @author rambaut
 *         Date: Apr 22, 2005
 *         Time: 9:12:27 AM
 */
public class Population {

    public Population(int populationSize, GenePool genePool) {

        this.populationSize = populationSize;
        this.genePool = genePool;
		this.selector = new NaturalSelector();

        lastGeneration = new Individual[populationSize];
        currentGeneration = new Individual[populationSize];

        Genome ancestor = genePool.createGenome(Genome.getMasterSequence());

        for (int i = 0; i < populationSize; i++) {
            currentGeneration[i] = new Individual(ancestor, null);
            lastGeneration[i] = new Individual();
        }

        ancestor.setFrequency(populationSize);

    }

    public void selectNextGeneration(Mutator mutator, boolean updateMaster) {

		double totalFitness = selector.initializeSelection(currentGeneration);
        meanFitness = totalFitness / populationSize;

        // first swap the arrays around
        Individual[] tmp = currentGeneration;
        currentGeneration = lastGeneration;
        lastGeneration = tmp;

        if (updateMaster) {
            genePool.resetStateFrequencies(populationSize);
        }

        // then select the currentGeneration based on the last.
        for (int i = 0; i < populationSize; i++) {
	        Individual parent = selector.selectIndividual();

            Genome genome = genePool.replicateGenome(parent.getGenome(), mutator);

            currentGeneration[i].setGenome(genome);
            currentGeneration[i].setParent(parent);

            if (updateMaster) {
                genePool.updateStateFrequencies(genome);
            }
        }

        // then kill off the virusexplorer in the last population.
        for (int i = 0; i < populationSize; i++) {
            genePool.killGenome(lastGeneration[i].getGenome());
        }

        if (updateMaster) {
            System.out.println("updating master");
            byte[] consensus = genePool.getConsensusSequence();
            genePool.changeMasterSequence(consensus);
        }
    }

    public void collectStatistics() {

        double d = 0;
        maxFrequency = 0;

        mostFrequentGenome = 0;

        for (int i = 0; i < populationSize; i++) {
            Genome genome = currentGeneration[i].getGenome();

            d += genome.getMutationCount();

            if (genome.getFrequency() > maxFrequency) {
                mostFrequentGenome = i;
                maxFrequency = genome.getFrequency();
            }
        }

        meanDistance = d / populationSize;
    }

    private final int populationSize;

    private final GenePool genePool;
	private final Selector selector;

    private Individual[] lastGeneration;
    private Individual[] currentGeneration;

    private int mostFrequentGenome;
    private int maxFrequency;
    private double meanDistance;

	private double meanFitness;

    public static void main(String[] args) {

        int genomeLength = 3300;
        int populationSize = 10000;
        int generationCount = 10;
        double mutationRate = 1E-4;

        Runtime rt = Runtime.getRuntime();

        double usedMemory = (rt.totalMemory() - rt.freeMemory()) / (1024*1024);
        System.out.println("Initial memory used: " + usedMemory + "MB");

        FitnessFunction aminoAcidFitnessFunction = new GammaFitnessFunction(genomeLength, 0.2, 1.0, 20, true);
        FitnessFunction fitness = new CodonFitnessFunction(genomeLength, GeneticCode.UNIVERSAL, aminoAcidFitnessFunction);

        byte[] masterSequence = fitness.getFittestSequence();

        Genome.setMasterSequence(masterSequence, fitness.getLogFitness(masterSequence));

        GenePool genePool = new GenePool(genomeLength, Codons.UNIVERSAL, fitness);

        System.out.print("Initializing population: ");
        Population population = new Population(populationSize, genePool);
        System.out.println("done.");

        usedMemory = (rt.totalMemory() - rt.freeMemory()) / (1024*1024);
        System.out.println("Base memory used: " + usedMemory + "MB");

        Mutator mutator = new CodonMutator(mutationRate, 2.0);

        for (int i = 0; i < generationCount; i++) {
            System.out.print("Generation " + i + ": ");

            boolean updateMaster = (i != 0 && i % 100 == 0);

            population.selectNextGeneration(mutator, updateMaster);

            population.collectStatistics();

            usedMemory = (rt.totalMemory() - rt.freeMemory()) / (1024*1024);

            System.out.println("fitness = " + population.meanFitness +
                                ", distance = " + population.meanDistance +
                                ", max freq = " + population.maxFrequency +
                                ", genepool size= " + genePool.getUniqueGenomeCount() + "(" + genePool.getUnusedGenomeCount() + " available) [" +
                                usedMemory + "MB]");
	        System.out.println(genePool.getDNASequenceString((Genome)genePool.getGenomes().get(population.mostFrequentGenome)));
	        System.out.println(genePool.getAminoAcidSequenceString((Genome)genePool.getGenomes().get(population.mostFrequentGenome)));
        }
    }

}
