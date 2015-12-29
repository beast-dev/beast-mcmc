/*
 * Population.java
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

import dr.math.MathUtils;

import java.util.ArrayList;
import java.util.List;

public class Population {

	// the sequences of this population
	List<Genome> population = null;

	double meanFitness;

	// mutates the sequences
	Mutator mutator = null;

	// calculates the fitness of the sequences
	FitnessFunction fitnessFunction = null;

	public Population(int size, int sequenceLength, Mutator mutator, FitnessFunction fitnessFunction, boolean initializeToFittest) {
		population = new ArrayList<Genome>();
		for (int i = 0; i < size; i++) {
			population.add(new SimpleGenome(sequenceLength, fitnessFunction, initializeToFittest));
		}
		this.mutator = mutator;
		this.fitnessFunction = fitnessFunction;
		cumulativeFitness = new double[size];
	}

	private Population(List<Genome> population, Mutator mutator, FitnessFunction fitnessFunction) {
		this.population = population;
		this.mutator = mutator;
		this.fitnessFunction = fitnessFunction;
	}

	/**
	 * pick parents randomly but proportional to fitness.
	 */
	public int pickParent() {

		return MathUtils.randomChoice(cumulativeFitness);
	}

	public final Mutator getMutator() { return mutator; }

	public final int getGenomeLength() {
		return population.get(0).getGenomeLength();
	}

	public final SimpleGenome getGenome(int i) {
		return (SimpleGenome)population.get(i);
	}

	/**
	 * This function returns the next generation!
	 */
	public final Population nextGeneration() {

		calculateCumulativeRelativeFitnesses();

		int popSize = population.size();

		List<Genome> nextPopulation = new ArrayList<Genome>();
		for (int i = 0; i < popSize; i++) {
			Genome parent = population.get(pickParent());
			nextPopulation.add(parent.replicate(mutator, fitnessFunction));
		}

		Population p = new Population(nextPopulation, mutator, fitnessFunction);

		return p;
	}

	//************************************************************************
	// private methods
	//************************************************************************

	private void calculateCumulativeRelativeFitnesses() {

		int popSize = population.size();

		if (cumulativeFitness == null) {
			cumulativeFitness = new double[popSize];
		}

		//System.err.println("popSize = " + popSize);
		//System.err.println("cumulativeFitness.length = " + cumulativeFitness.length);


		cumulativeFitness[0] =  ((SimpleGenome)population.get(0)).fitness;

		for (int i = 1; i < popSize; i++) {
			cumulativeFitness[i] = cumulativeFitness[i-1] + ((SimpleGenome)population.get(i)).fitness;
		}
		double totalFitness = cumulativeFitness[popSize-1];
		if (totalFitness == 0.0) throw new RuntimeException("Population crashed! No viable children.");

		for (int i = 0; i < popSize; i++) {
			cumulativeFitness[i] /= totalFitness;
		}
		meanFitness = totalFitness / popSize;
	}


	//************************************************************************
	// static methods
	//************************************************************************


	public static Population forwardSimulation(Population startingPopulation, int generations) {
		Population p = startingPopulation;
		for (int i = 0; i < generations; i++) {
			p = p.nextGeneration();
		}
		return p;
	}

	public int getAgeOfMRCA(double[] f) {

		int popSize = population.size();

		for (int i = 0; i < popSize; i++) {
			getGenome(i).mark();
		}

		int age = 1;
		Genome g = getGenome(0);
		while (g != null && g.getMarks() < popSize) {
			g = g.getParent();
			age += 1;
		}

		for (int i = 0; i < popSize; i++) {
			getGenome(i).unmark();
		}

		if (g != null) {
			f[0] = g.getFitness();
			return age;
		} else {
			f[0] = 0;
			return -1;
		}
	}

	public Genome getMRCA() {

		int popSize = population.size();

		for (int i = 0; i < popSize; i++) {
			getGenome(i).mark();
		}

		Genome g = getGenome(0);
		while (g != null && g.getMarks() < popSize) {
			g = g.getParent();
		}

		for (int i = 0; i < popSize; i++) {
			getGenome(i).unmark();
		}

		return g;
	}

	public void unfoldedSiteFrequencies(List<Double>[] siteFrequencies) {

		final int popSize = population.size();

		double[][] fitnessTable = fitnessFunction.getFitnessTable();

		Genome progenitorGenome = getMRCA();
		byte[] progenitor = progenitorGenome.getSequence();

		int genomeLength = getGenome(0).getGenomeLength();

		for (int j = 0; j < genomeLength; j++) {

			byte wildType = progenitor[j];
			double wildFitness = fitnessTable[j][wildType];
			double mutantFitness = fitnessTable[j][(wildType+1)%2];
			double relativeFitness = mutantFitness/wildFitness;

			if (wildFitness == 0.0) {
				for (int k = 0; k < genomeLength; k++) {
					System.err.print(progenitor[k]);
				}
				System.err.println();
				throw new RuntimeException("Progenitor fitness can not be zero!");
			}

			if (relativeFitness != 0.0) {
				int mutantCount = 0;
				for (int i = 0; i < popSize; i++) {
					if (getGenome(i).sequence[j] != wildType) {
						mutantCount += 1;
					}
				}
				siteFrequencies[mutantCount].add(relativeFitness);
			}
		}
	}

	public double getMeanParentFitness() {
		double parentFitness = 0.0;
		int popSize = population.size();
		for (int i = 0; i < popSize; i++) {
			parentFitness += getGenome(i).getParent().getFitness();
		}
		return parentFitness/(double)popSize;
	}

	public double getProportionAsFit(double fit) {
		int count = 0;
		int popSize = population.size();
		for (int i = 0; i < popSize; i++) {
			if (getGenome(i).fitness >= fit) count += 1;
		}
		return (double)count/(double)popSize;
	}

	public void getTipHammingDistance(int individuals, List<Double>[] distances) {

		int popSize = population.size();

		for (int i = 0; i < individuals; i++) {
			getGenome(i).mark();
		}

		for (int i = 0; i < individuals; i++) {
			Genome tip = getGenome(i);
			Genome a = tip;
			int generation = 0;
			while ((a.getMarks() < popSize) && (generation < distances.length)) {
				distances[i].add((double) tip.hammingDistance(a));
				generation += 1;
				a = a.getParent();
			}
		}

		for (int i = 0; i < individuals; i++) {
			getGenome(i).unmark();
		}
	}

	public void getMutationDensity(int individuals, List<Integer>[] mutations) {

		for (int i = 0; i < individuals; i++) {
			getGenome(i).mark();
		}

		for (int i = 0; i < individuals; i++) {
			Genome genome = getGenome(i);
			int generation = 0;
			while ((genome.getMarks() < individuals) && (generation < mutations.length)) {
				mutations[generation].add(genome.getMutations().length);
				generation += 1;
				genome = genome.getParent();
			}
		}

		for (int i = 0; i < individuals; i++) {
			getGenome(i).unmark();
		}
	}

	public int getFirstActiveLineage() {
		int popSize = population.size();
		for (int i = 0; i < popSize; i++) {
			if (getGenome(i).getMarks() > 0) return i;
		}
		return -1;
	}


	public int getActiveLineageCount() {
		int count = 0;
		int popSize = population.size();
		for (int i = 0; i < popSize; i++) {
			if (getGenome(i).getMarks() > 0) count += 1;
		}
		return count;
	}

	public int getActiveMutationsCount() {
		int count = 0;
		int popSize = population.size();
		for (int i = 0; i < popSize; i++) {
			if (getGenome(i).getMarks() > 0) count += getGenome(i).mutations.length;
		}
		return count;
	}



	static double[] cumulativeFitness = null;

}