/*
 * Mutator.java
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

import dr.math.matrixAlgebra.Matrix;

/**
 * @author Alexei Drummond
 */
public abstract class Mutator {

	Mutator() {}

	public Mutator(Matrix mutationRates) {
		this.mutationRates = mutationRates;
	}

	public Mutator(double mutationRate, int stateSize) {
		double noMutation = 1.0 - ((stateSize-1)*mutationRate);

		double[][] rates = new double[stateSize][stateSize];
		for (int i = 0; i < stateSize; i++) {
			for (int j = 0; j < stateSize; j++) {
				if (i != j) {
					rates[i][j] = mutationRate;
				} else rates[i][i] = noMutation;
			}
		}
		this.mutationRates = new Matrix( rates);
	}

	public abstract Mutation[] mutate(byte[] sequence, byte[] childSequence);

	Matrix mutationRates = null;
}