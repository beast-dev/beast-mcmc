/*
 * MultiChainRunner.java
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

package dr.inference.mcmcmc;

import java.util.List;
import java.util.concurrent.*;
import java.util.concurrent.atomic.AtomicInteger;

/**
 * Runs one or more chain sets to completion. Each set is driven from its own thread and the
 * sets never wait for each other. Every chain has its own thread, and at most chainThreads
 * chains compute at any one time.
 *
 * If any set fails, all the others are asked to stop and the first failure is rethrown.
 *
 * @author Andrew Rambaut
 */
public class MultiChainRunner {

    public MultiChainRunner(List<ChainSet> sets, int chainThreads) {
        if (sets.isEmpty()) {
            throw new IllegalArgumentException("No chains to run");
        }
        if (chainThreads < 1) {
            throw new IllegalArgumentException("At least one chain thread is needed");
        }
        this.sets = sets;
        this.chainThreads = chainThreads;
    }

    public void run() {
        final Semaphore permits = new Semaphore(chainThreads, true);

        final AtomicInteger driverCount = new AtomicInteger(0);
        ExecutorService drivers = Executors.newFixedThreadPool(sets.size(),
                r -> new Thread(r, "chain-set-" + driverCount.incrementAndGet()));
        CompletionService<ChainSet> completion = new ExecutorCompletionService<>(drivers);

        timer.start();

        for (final ChainSet set : sets) {
            completion.submit(() -> {
                set.start(permits);
                while (set.step()) {
                    // the set reports its own progress through its loggers
                }
                set.finish();
                return set;
            });
        }

        Throwable failure = null;
        try {
            for (int i = 0; i < sets.size(); i++) {
                try {
                    completion.take().get();
                } catch (ExecutionException ee) {
                    if (failure == null) {
                        failure = ee.getCause();
                        pleaseStop();
                    }
                } catch (InterruptedException ie) {
                    // interrupted from outside (e.g., the console window): stop and wait for
                    // the sets to wind down
                    pleaseStop();
                    i--;
                }
            }
        } finally {
            for (ChainSet set : sets) {
                set.shutdown();
            }
            drivers.shutdownNow();
            timer.stop();
        }

        System.out.println();
        System.out.println("Time taken: " + timer.toString());

        if (failure instanceof RuntimeException) {
            throw (RuntimeException) failure;
        } else if (failure instanceof Error) {
            throw (Error) failure;
        } else if (failure != null) {
            throw new RuntimeException(failure);
        }
    }

    /**
     * Requests that every chain stops. Can be called from any thread.
     */
    public void pleaseStop() {
        for (ChainSet set : sets) {
            set.pleaseStop();
        }
    }

    public List<ChainSet> getChainSets() {
        return sets;
    }

    public dr.util.Timer getTimer() {
        return timer;
    }

    private final List<ChainSet> sets;
    private final int chainThreads;
    private final dr.util.Timer timer = new dr.util.Timer();
}
