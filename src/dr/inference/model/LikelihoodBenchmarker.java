/*
 * LikelihoodBenchmarker.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
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
 *
 */

package dr.inference.model;

import java.util.Arrays;
import java.util.List;
import java.util.logging.Logger;

/**
 * @author Andrew Rambaut
 * @author Frederik M. Andersen
 */
public class LikelihoodBenchmarker {

    public LikelihoodBenchmarker(List<Likelihood> likelihoods, int iterationCount) {
        Logger logger = Logger.getLogger("dr.app.beagle");

        for (Likelihood likelihood : likelihoods) {
            double[] times = new double[iterationCount];

            for (int i = 0; i < iterationCount; i++) {
                likelihood.makeDirty();
                long startTime = System.nanoTime();
                likelihood.getLogLikelihood();
                long endTime = System.nanoTime();
                times[i] = (endTime - startTime) * 1E-3; // microseconds
            }

            Arrays.sort(times);

            double sum = 0.0;
            for (double t : times) {
                sum += t;
            }
            double mean = sum / iterationCount;

            double sumSqDev = 0.0;
            for (double t : times) {
                double dev = t - mean;
                sumSqDev += dev * dev;
            }
            double sd = Math.sqrt(sumSqDev / iterationCount);

            double median = iterationCount % 2 == 0
                    ? (times[iterationCount / 2 - 1] + times[iterationCount / 2]) / 2.0
                    : times[iterationCount / 2];

            double min = times[0];
            double max = times[iterationCount - 1];

            String name = likelihood.getId() + " (" + likelihood.getClass().getSimpleName() + ")";
            logger.info(String.format("Benchmark %s [%d iterations]: " +
                            "mean=%.1f us, sd=%.1f us, median=%.1f us, range=[%.1f, %.1f] us",
                    name, iterationCount, mean, sd, median, min, max));
        }
    }
}
