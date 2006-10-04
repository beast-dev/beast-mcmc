/*
 * GammaFitnessFunction.java
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

package dr.simulator.fitnessfunctions;

import dr.math.GammaDistribution;
import dr.simulator.random.Random;

public class GammaFitnessFunction extends AbstractFitnessFunction {

	/**
	 * selection coefficients are gamma distributed with parameters alpha and beta
	 */
	public GammaFitnessFunction(int genomeLength, double alpha, double beta, int stateSize, boolean randomFittest) {

        double[][] logFitnessTable = new double[genomeLength][stateSize];

		GammaDistribution gamma = new GammaDistribution(beta, alpha);

		int fitpos = 0;
		for (int i = 0; i < genomeLength; i++) {

			if (randomFittest) {
				fitpos = Random.nextInt(stateSize);
			}
			logFitnessTable[i][fitpos] = 0.0;
			for (int j = 0; j < stateSize; j++) {

				if (j != fitpos) {
					double prob = Math.round(Random.nextDouble() * 1000.0)/1000.0;
					while ((prob <= 0.0) || (prob >= 1.0)) {
						prob = Math.round(Random.nextDouble() * 1000.0)/1000.0;
					}
					logFitnessTable[i][j] = Math.log(Math.max(0.0, 1.0 - gamma.quantile(prob)));
				}
			}
		}

        setLogFitnessTable(logFitnessTable);
	}

}