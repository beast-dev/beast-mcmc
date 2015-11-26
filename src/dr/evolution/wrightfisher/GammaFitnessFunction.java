/*
 * GammaFitnessFunction.java
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
import dr.math.distributions.GammaDistribution;

public class GammaFitnessFunction extends FitnessFunction {

    /**
     * selection coefficients are gamma distributed with parameters alpha and beta
     */
    public GammaFitnessFunction(int genomeLength, double alpha, double beta, int stateSize, boolean randomFittest) {

        GammaDistribution gamma = new GammaDistribution(beta, alpha);

        fitness = new double[genomeLength][stateSize];
        fittest = new byte[genomeLength];
        int fitpos = 0;
        for (int i = 0; i < genomeLength; i++) {

            if (randomFittest) {
                fitpos = MathUtils.nextInt(stateSize);
            }
            fitness[i][fitpos] = 1.0;
            fittest[i] = (byte) fitpos;
            for (int j = 0; j < stateSize; j++) {

                if (j != fitpos) {
                    double prob = Math.round(MathUtils.nextDouble() * 1000.0) / 1000.0;
                    while ((prob <= 0.0) || (prob >= 1.0)) {
                        prob = Math.round(MathUtils.nextDouble() * 1000.0) / 1000.0;
                    }
                    fitness[i][j] = Math.max(0.0, 1.0 - gamma.quantile(prob));
                }
            }
        }
    }

    public final double getFitness(byte[] sequence) {

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

        for (int j = 0; j < fitness[0].length; j++) {
            for (int i = 0; i < fitness.length; i++) {
                System.out.print((Math.round(fitness[i][j] * 1000.0) / 1000.0) + "\t");
            }
            System.out.println();
        }

        return fitness;
    }

    public void initializeToFittest(byte[] genome) {
        System.arraycopy(fittest, 0, genome, 0, fittest.length);
    }

    double[][] fitness;
    byte[] fittest = null;
}