/*
 * CodonMutator.java
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
public class CodonMutator implements Mutator {

	public CodonMutator(double muRate, double transitionTransversionRatio) {
        this.transitionProbability = transitionTransversionRatio / (1.0 + transitionTransversionRatio);
		this.muRate = muRate;

        codonTable = new byte[64][];
        int u = 0;
        for (byte i = 0; i < 4; i++) {
            for (byte j = 0; j < 4; j++) {
                for (byte k = 0; k < 4; k++) {
                    codonTable[u] = new byte[] { i, j, k };
                    u++;
                }
            }
        }
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

            byte oldState = genome.getState(pos);

            byte[] codon = codonTable[oldState];

            byte newState;

            int codonPosition = Random.nextInt(3);
            switch (codonPosition) {
                case 0:
                    newState = (byte)((mutateNucleotide(codon[0]) * 16) + (codon[1] * 4) + codon[2]);
                break;
                case 1:
                    newState = (byte)((codon[0] * 16) + (mutateNucleotide(codon[1]) * 4) + codon[2]);
                break;
                case 2:
                    newState = (byte)((codon[0] * 16) + (codon[1] * 4) + mutateNucleotide(codon[2]));
                break;
                default: newState = -1;
            }

			mutations[i] = Mutation.getMutation(pos, newState);
		}

        return mutations;
	}

    byte mutateNucleotide(byte oldState) {
        boolean isTransition = Random.nextBoolean(transitionProbability);
        if (isTransition) {
            switch (oldState) {
                case Nucleotides.A_STATE: return Nucleotides.G_STATE;
                case Nucleotides.C_STATE: return Nucleotides.UT_STATE;
                case Nucleotides.G_STATE: return Nucleotides.A_STATE;
                case Nucleotides.UT_STATE: return Nucleotides.C_STATE;
            }
        } else {
            switch (oldState) {
                case Nucleotides.A_STATE:
                case Nucleotides.G_STATE:
                    return (byte)(Random.nextBoolean() ? Nucleotides.C_STATE : Nucleotides.UT_STATE);
                case Nucleotides.C_STATE:
                case Nucleotides.UT_STATE:
                    return (byte)(Random.nextBoolean() ? Nucleotides.A_STATE : Nucleotides.G_STATE);
            }
        }
        return -1;
    }

	private final double muRate;
    private final double transitionProbability;
	private final byte[][] codonTable;

	private double poissonMean = -1;
	private int genomeLength = -1;
}