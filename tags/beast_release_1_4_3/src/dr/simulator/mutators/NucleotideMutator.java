/*
 * NucleotideMutator.java
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

package dr.simulator.mutators;

import dr.simulator.Mutation;
import dr.simulator.Genome;
import dr.simulator.datatypes.Nucleotides;
import dr.simulator.random.Random;

/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 */
public class NucleotideMutator implements Mutator {

	public NucleotideMutator(double muRate, double transitionTransversionRatio) {
        this.transitionProbability = transitionTransversionRatio / (1.0 + transitionTransversionRatio);
		this.muRate = muRate;
	}

	public Mutation[] mutate(Genome genome) {

		if (genomeLength != genome.getLength()) {
			genomeLength = genome.getLength();
			poissonMean = genomeLength * muRate;
		}

		int mutationCount = Random.nextPoisson(poissonMean);

        Mutation[] mutations = new Mutation[mutationCount];

		for (int i = 0; i < mutationCount; i++) {
			int pos = Random.nextInt(genomeLength);
            boolean isTransition = Random.nextBoolean(transitionProbability);
            byte oldState = genome.getState(pos);
            byte newState;

            if (isTransition) {
                switch (oldState) {
                    case Nucleotides.A_STATE: newState = Nucleotides.G_STATE; break;
                    case Nucleotides.C_STATE: newState = Nucleotides.UT_STATE; break;
                    case Nucleotides.G_STATE: newState = Nucleotides.A_STATE; break;
                    case Nucleotides.UT_STATE: newState = Nucleotides.C_STATE; break;
                    default: newState = -1; // should never get this
                }
            } else {
                switch (oldState) {
                    case Nucleotides.A_STATE:
                    case Nucleotides.G_STATE:
                        newState = (byte)(Random.nextBoolean() ? Nucleotides.C_STATE : Nucleotides.UT_STATE);
                        break;
                    case Nucleotides.C_STATE:
                    case Nucleotides.UT_STATE:
                        newState = (byte)(Random.nextBoolean() ? Nucleotides.A_STATE : Nucleotides.G_STATE);
                        break;
                    default: newState = -1; // should never get this
                }
            }

			mutations[i] = Mutation.getMutation(pos, newState);
		}

        return mutations;
	}

	private final double muRate;
    private final double transitionProbability;

	private double poissonMean = -1;
	private int genomeLength = -1;
}