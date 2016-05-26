/*
 * BakSneppenFitness.java
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

public class BakSneppenFitness extends FitnessFunction {
	
	public BakSneppenFitness(double M) {
		if (M >= 1.0) this.M = M;
		else throw new IllegalArgumentException("M must be greater than or equal to 1");
	}
	
	public double getFitness(byte[] sequence) {
		double fitness = 0.0;
		double currentFitness = 1.0;
		for (int i = 0; i < sequence.length; i++) {
			if (sequence[i] == 0) fitness += currentFitness;
			currentFitness /= M;
		}
		return fitness;
	}
	
	public double getFitnessFactor(int pos, byte state1, byte state2) {
		throw new UnsupportedOperationException();
	}
	
	public double[][] getFitnessTable() {
		throw new UnsupportedOperationException();
	}
	
	public void initializeToFittest(byte[] genome) {
		throw new UnsupportedOperationException();
	}
	
	private double M = 2.0;
}