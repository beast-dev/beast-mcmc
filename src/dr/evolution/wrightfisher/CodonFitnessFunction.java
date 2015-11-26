/*
 * CodonFitnessFunction.java
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

import dr.evolution.datatype.AminoAcids;
import dr.evolution.datatype.Codons;
import dr.evolution.datatype.GeneticCode;
import dr.math.MathUtils;
import dr.math.distributions.GammaDistribution;

public class CodonFitnessFunction extends FitnessFunction {

    /**
     * selection coefficients are gamma distributed with parameters alpha and beta
     */
    public CodonFitnessFunction(int codonCount, double alpha, double beta, double pInv) {

        GammaDistribution gamma = new GammaDistribution(beta, alpha);
        stateSize = 20;

        fitness = new double[codonCount][stateSize];
        byte[] aaFittest = new byte[codonCount];
        fittest = getFittestGenome(codonCount, aaFittest);
        for (int i = 0; i < codonCount; i++) {

            fitness[i][aaFittest[i]] = 1.0;
            for (int j = 0; j < stateSize; j++) {

                if (j != aaFittest[i]) {

                    if (MathUtils.nextDouble() < pInv) {
                        fitness[i][j] = 0.0;
                    } else {

                        double prob = Math.round(MathUtils.nextDouble() * 1000.0) / 1000.0;
                        while ((prob <= 0.0) || (prob >= 1.0)) {
                            prob = Math.round(MathUtils.nextDouble() * 1000.0) / 1000.0;
                        }
                        fitness[i][j] = Math.max(0.0, 1.0 - gamma.quantile(prob));
                    }
                }
            }
        }
    }

    /**
     * Assumes the sequence is represented as nucleotides
     */
    public final double getFitness(byte[] sequence) {

        double totalFitness = 1.0;

        int count = 0;
        for (int i = 0; i < sequence.length; i += 3) {
            int codonState = (sequence[i] * 16) + (sequence[i + 1] * 4) + sequence[i + 2];
            //int codonState = codons.getState(sequence[i], sequence[i+1], sequence[i+2]);

            int aminoAcid = geneticCode.getAminoAcidState(codonState);

            //System.out.println(sequence[i] + "" + sequence[i+1] + "" + sequence[i+2] +"="+ AminoAcids.INSTANCE.getChar(aminoAcid));

            // stop codons are lethal of course
            if (aminoAcid < 0 || aminoAcid >= stateSize) {
                //System.out.println("Stop codon found! " + aminoAcid);
                return 0.0;
            }
            totalFitness *= fitness[count][aminoAcid];
            if (totalFitness == 0.0) return 0.0;
            count += 1;
        }
        return totalFitness;
    }

    /**
     * @return the relative fitness increase of the new state at given position to the old state.
     */
    public double getFitnessFactor(int pos, byte newState, byte oldState) {
        throw new RuntimeException();
    }

    public final double[][] getFitnessTable() {

        for (int j = 0; j < fitness[0].length; j++) {
            for (int i = 0; i < fitness.length; i++) {
                System.out.print((Math.round(fitness[i][j] * 1000.0) / 1000.0) + "\t");
            }
            System.out.println();
        }

        return fitness;
    }

    public byte[] getFittestGenome(int codonCount, byte[] aa) {

        byte[] genome = new byte[codonCount * 3];

        for (int i = 0; i < codonCount; i++) {

            int aminoAcid = stateSize;
            int pos1 = MathUtils.nextInt(4);
            int pos2 = MathUtils.nextInt(4);
            int pos3 = MathUtils.nextInt(4);
            int codonState = (pos1 * 16) + (pos2 * 4) + pos3;
            //int codonState = codons.getState(pos1, pos2, pos3);
            aminoAcid = geneticCode.getAminoAcidState(codonState);

            while (aminoAcid >= stateSize) {
                pos1 = MathUtils.nextInt(4);
                pos2 = MathUtils.nextInt(4);
                pos3 = MathUtils.nextInt(4);
                codonState = (pos1 * 16) + (pos2 * 4) + pos3;
                //codonState = codons.getState(pos1, pos2, pos3);
                aminoAcid = geneticCode.getAminoAcidState(codonState);
            }
            System.out.print(AminoAcids.INSTANCE.getChar(aminoAcid));

            genome[i * 3] = (byte) pos1;
            genome[i * 3 + 1] = (byte) pos2;
            genome[i * 3 + 2] = (byte) pos3;
            aa[i] = (byte) aminoAcid;
        }
        System.out.println();

        return genome;
    }

    public void initializeToFittest(byte[] genome) {
        System.arraycopy(fittest, 0, genome, 0, fittest.length);
    }

    int stateSize = 20;
    double[][] fitness;
    byte[] fittest = null;
    GeneticCode geneticCode = GeneticCode.UNIVERSAL;
    Codons codons = Codons.UNIVERSAL;
}