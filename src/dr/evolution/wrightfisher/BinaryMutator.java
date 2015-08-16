/*
 * BinaryMutator.java
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
import dr.math.Poisson;

/**
 * @author Alexei Drummond
 */
public class BinaryMutator extends Mutator {

	public BinaryMutator(double muRate) {

		this.muRate = muRate;
	}

	public Mutation[] mutate(byte[] sequence, byte[] childSequence) {

		if (genomeLength != sequence.length) {
			genomeLength = sequence.length;
			poissonMean = genomeLength * muRate;
		}

		int mutationCount = Poisson.nextPoisson(poissonMean);
		Mutation[] mutations = new Mutation[mutationCount];

		System.arraycopy(sequence, 0, childSequence, 0, genomeLength);

		for (int i = 0; i < mutationCount; i++) {
			int pos = MathUtils.nextInt(genomeLength);
			byte state = (byte)((sequence[pos] + 1) % 2);
			childSequence[pos] = state;
			mutations[i] = new Mutation(pos, state);
		}
		return mutations;
	}

	double muRate = 0.01;
	double poissonMean = -1;
	int genomeLength = -1;
}