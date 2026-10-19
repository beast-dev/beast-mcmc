/*
 * AbstractChainSet.java
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

import dr.inference.mcmc.ChainContext;
import dr.inference.mcmc.MCMC;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.*;

/**
 * Owns one thread per chain so that each chain (and its BEAGLE instances) is always run
 * from the same thread, with the chain's context (and so its random number stream) bound to
 * it. A chain computes only while it holds a permit from the semaphore shared across the run.
 *
 * @author Andrew Rambaut
 */
public abstract class AbstractChainSet implements ChainSet {

    /**
     * @param contexts      one per chain, bound to the chain's thread
     * @param driverContext bound to the thread driving the set (e.g., for swap decisions)
     */
    protected AbstractChainSet(int index, MCMC[] mcmcs, ChainContext[] contexts, ChainContext driverContext) {
        if (contexts.length != mcmcs.length) {
            throw new IllegalArgumentException("There should be one context per chain");
        }
        this.index = index;
        this.mcmcs = mcmcs;
        this.contexts = contexts;
        this.driverContext = driverContext;
    }

    @Override
    public int getIndex() {
        return index;
    }

    @Override
    public void start(Semaphore permits) {
        // start is called on the driver thread
        driverContext.bind();

        this.permits = permits;
        workers = new ExecutorService[mcmcs.length];
        for (int i = 0; i < mcmcs.length; i++) {
            final String name = "chain-" + (index + 1) + (mcmcs.length > 1 ? "." + i : "");
            final ChainContext context = contexts[i];
            workers[i] = Executors.newSingleThreadExecutor(r -> new Thread(() -> {
                context.bind();
                r.run();
            }, name));
        }
    }

    @Override
    public void shutdown() {
        if (workers != null) {
            for (ExecutorService worker : workers) {
                worker.shutdownNow();
            }
        }
    }

    @Override
    public void pleaseStop() {
        for (MCMC mcmc : mcmcs) {
            mcmc.pleaseStop();
        }
    }

    /**
     * Runs one task per chain, each on that chain's own thread, and waits for all of them.
     * Any exception from a chain is rethrown here after the others have finished.
     */
    protected <T> List<T> runOnChains(List<? extends Callable<T>> tasks) {
        List<Future<T>> futures = new ArrayList<>();
        for (int i = 0; i < tasks.size(); i++) {
            final Callable<T> task = tasks.get(i);
            futures.add(workers[i].submit(() -> {
                permits.acquire();
                try {
                    return task.call();
                } finally {
                    permits.release();
                }
            }));
        }

        List<T> results = new ArrayList<>();
        Throwable failure = null;
        for (Future<T> future : futures) {
            try {
                results.add(getUninterruptibly(future));
            } catch (ExecutionException ee) {
                if (failure == null) {
                    failure = ee.getCause();
                    // stop the rest of the set promptly
                    pleaseStop();
                }
            }
        }

        if (failure instanceof RuntimeException) {
            throw (RuntimeException) failure;
        } else if (failure instanceof Error) {
            throw (Error) failure;
        } else if (failure != null) {
            throw new RuntimeException(failure);
        }
        return results;
    }

    private static <T> T getUninterruptibly(Future<T> future) throws ExecutionException {
        boolean interrupted = false;
        try {
            while (true) {
                try {
                    return future.get();
                } catch (InterruptedException ie) {
                    interrupted = true;
                }
            }
        } finally {
            if (interrupted) {
                Thread.currentThread().interrupt();
            }
        }
    }

    protected final int index;
    protected final MCMC[] mcmcs;
    protected final ChainContext[] contexts;
    private final ChainContext driverContext;

    private Semaphore permits;
    private ExecutorService[] workers;
}
