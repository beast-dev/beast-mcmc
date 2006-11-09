/*
 * NaturalSelector.java
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

package dr.simulator.selectors;

import dr.simulator.Individual;
import dr.simulator.random.Random;

import java.util.Arrays;

/**
 * @author rambaut
 *         Date: Apr 27, 2005
 *         Time: 2:50:36 PM
 */
public class NaturalSelector implements Selector {

    public double initializeSelection(Individual[] population) {

        this.population = population;

        if (populationSize != population.length) {
            populationSize = population.length;
            cumulativeFitness = new double[populationSize];
	        randomNumbers = new double[populationSize];
	    }

        cumulativeFitness[0] = Math.exp(population[0].getLogFitness());

	    for (int i = 1; i < populationSize; i++) {
		    cumulativeFitness[i] = cumulativeFitness[i-1] + Math.exp(population[i].getLogFitness());
	    }

	    double totalFitness = cumulativeFitness[populationSize-1];
	    if (totalFitness == 0.0) throw new RuntimeException("Population crashed! No viable children.");

	    for (int i = 0; i < populationSize; i++) {
		    cumulativeFitness[i] /= totalFitness;
	    }

	    for (int i = 0; i < populationSize; i++) {
		    randomNumbers[i] = Random.nextDouble();
	    }

        Arrays.sort(randomNumbers);

        currentIndividual = 0;
        currentRandomNumber = 0;

        return totalFitness;
    }

    /**
     * pick parents randomly but proportional to fitness.
     */
    public Individual selectIndividual() {

        if (currentIndividual >= populationSize || currentRandomNumber >= populationSize) {
            throw new IllegalArgumentException("No more individuals available to select");
        }

        while (randomNumbers[currentRandomNumber] > cumulativeFitness[currentIndividual]) {
			currentIndividual++;
        }
        currentRandomNumber++;

	    return population[currentIndividual];
    }

    private Individual[] population;
    private int populationSize = 0;
    private int currentIndividual = -1;
    private int currentRandomNumber = -1;
    private double[] randomNumbers = null;
	private double[] cumulativeFitness = null;
}
