/*
 * ParallelRobustCountingExecutor.java
 *
 * Copyright (c) 2002-2023 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodel.substmodel;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.*;

/**
 * @author Marc A. Suchard
 */

public class ParallelRobustCountingExecutor {

    private final ExecutorService pool;
    private final List<DerivativeCaller> derivativeCaller = new ArrayList<>();

    public ParallelRobustCountingExecutor(int requestedThreads, MarkovJumpsSubstitutionModel model,
                                          int numSites) {

        int threadCount = Math.min(requestedThreads, numSites);
        pool = Executors.newFixedThreadPool(threadCount);

        int elementsPerThread = (numSites + threadCount - 1) / threadCount;
        int start = 0;
        for (int i = 0; i < threadCount; ++i) {
            derivativeCaller.add(new DerivativeCaller(start, Math.min(start + elementsPerThread, numSites), model));
            start += threadCount;
        }
    }

    public void execute(double[] out) {

        for (DerivativeCaller caller : derivativeCaller) {
            caller.setOutput(out);
        }

        try {
            pool.invokeAll(derivativeCaller);
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
    }

    private static class DerivativeCaller implements Callable<Void> {

        private final int start;
        private final int end;

        private double[] out;

        public DerivativeCaller(int start, int end, MarkovJumpsSubstitutionModel model) {
            this.start = start;
            this.end = end;
        }

        public void setOutput(double[] out) {
            this.out = out;
        }

        public Void call() throws Exception {
            for (int i = start; i < end; ++i) {
                out[i] = 0.0;
            }
            return null;
        }
    }

    public static final boolean DEBUG_PARALLEL_EVALUATION = false;
}
