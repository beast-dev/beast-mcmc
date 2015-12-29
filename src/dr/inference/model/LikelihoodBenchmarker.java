/*
 * LikelihoodBenchmarker.java
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

package dr.inference.model;

import java.util.List;
import java.util.logging.Logger;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public class LikelihoodBenchmarker {

    public LikelihoodBenchmarker(List<Likelihood> likelihoods, int iterationCount) {
        for (Likelihood likelihood : likelihoods) {
            long startTime = System.nanoTime();

            for (int i = 0; i < iterationCount; i++) {
                likelihood.makeDirty();
                likelihood.getLogLikelihood();
            }

            long endTime = System.nanoTime();

            double seconds = (endTime - startTime) * 1E-9;
            Logger.getLogger("dr.app.beagle").info(
                    "Benchmark " + likelihood.getId() + "(" + likelihood.getClass().getName() + "): " +
                            seconds + " sec");

        }
    }
}
