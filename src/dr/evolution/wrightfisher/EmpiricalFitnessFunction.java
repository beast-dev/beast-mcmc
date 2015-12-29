/*
 * EmpiricalFitnessFunction.java
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

public class EmpiricalFitnessFunction extends FitnessFunction {

	public EmpiricalFitnessFunction(int genomeLength, double[] mutationFitnesses, int stateSize, boolean isRandom) {
		
		fitness = new double[genomeLength][stateSize];
		logFitness = new double[genomeLength][stateSize];
		int index = 0;
		for (int i = 0; i < genomeLength; i++) {
			fitness[i][0] = 1.0;
			for (int j = 1; j < stateSize; j++) {
				if (isRandom) {
					fitness[i][j] = mutationFitnesses[MathUtils.nextInt(mutationFitnesses.length)];
				} else {
					fitness[i][j] = mutationFitnesses[index];
					index = (index + 1) % mutationFitnesses.length;
				}
			}
		}
		
		// tabulate log fitnesses
		for (int i = 0; i < stateSize; i++) {
			for (int j = 0; j < genomeLength; j++) {
				logFitness[j][i] = Math.log(fitness[j][i]);
			}
		}
	}
	
	public final double getFitness(byte[] sequence) {
		
		/*double totalFitness = 0.0;
		
		for (int i = 0; i < sequence.length; i++) {
			totalFitness += logFitness[i][sequence[i]];
		}
		return Math.exp(totalFitness);
		*/
		double totalFitness = 1.0;
		
		for (int i = 0; i < sequence.length; i++) {
			totalFitness *= fitness[i][sequence[i]];
		}
		return totalFitness;
	}
	
	/**
	 * @return the relative fitness increase of the new state at given position to the old state.
	 */
	public double getFitnessFactor(int pos, byte newState, byte oldState) {
		return fitness[pos][newState] / fitness[pos][oldState];
	}
	
	public final double[][] getFitnessTable() {
		return fitness;
	}
	
	public void initializeToFittest(byte[] genome) {
		for (int i = 0; i < genome.length; i++) {
			genome[i] = 0;
		}
	}

	
	double[][] fitness;
	double[][] logFitness;
}