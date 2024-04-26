/*
 * ParallelGradientExecutor.java
 *
 * Copyright (c) 2002-2022 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.inference.hmc;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.*;

/**
 * @author Marc A. Suchard
 */

public class ParallelGradientExecutor {

    private final ExecutorService pool;
    private final List<DerivativeCaller> derivativeCaller;

    interface Reducer {
        double[] reduce(List<Future<double[]>> gradients, int length) throws ExecutionException, InterruptedException;
    }

    public ParallelGradientExecutor(int threads, List<GradientWrtParameterProvider> derivativeList) {

        assert derivativeList.size() > 1;

        if (threads <= 0) {
            pool = Executors.newCachedThreadPool();
        } else {
            int threadCount = Math.min(threads, derivativeList.size());
            pool = Executors.newFixedThreadPool(threadCount);
        }

        derivativeCaller = new ArrayList<>(derivativeList.size());
        for (int i = 0; i < derivativeList.size(); ++i) {
            derivativeCaller.add(new DerivativeCaller(derivativeList.get(i), i));
        }
    }

    public double[] getDerivativeLogDensityInParallel(JointGradient.DerivativeType derivativeType,
                                                      Reducer reducer, int length) {

        for (DerivativeCaller caller : derivativeCaller) {
            caller.setDerivativeType(derivativeType);
        }

        double[] derivative = null;

        try {
            List<Future<double[]>> results = pool.invokeAll(derivativeCaller);
            derivative = reducer.reduce(results, length);
        } catch (InterruptedException | ExecutionException e) {
            e.printStackTrace();
        }

        return derivative;
    }

    private static class DerivativeCaller implements Callable<double[]> {

        public DerivativeCaller(GradientWrtParameterProvider gradient, int index) {
            this.gradient = gradient;
            this.index = index;
        }

        public double[] call() throws Exception {
            if (DEBUG_PARALLEL_EVALUATION) {
                System.err.println("Invoking thread #" + index + " for " + gradient.getLikelihood().getId() +
                        " with type " + type + " in executor");
            }

            return type.getDerivativeLogDensity(gradient);
        }

        public void setDerivativeType(JointGradient.DerivativeType type) {
            this.type = type;
        }

        private final GradientWrtParameterProvider gradient;
        private final int index;

        private JointGradient.DerivativeType type;
    }

    public static final boolean DEBUG_PARALLEL_EVALUATION = false;
}
