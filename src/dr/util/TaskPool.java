/*
 * TaskPool.java
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

package dr.util;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.*;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.locks.LockSupport;
import java.util.function.BinaryOperator;

/**
 * @author Marc A. Suchard
 * @author Andrew Holbrook
 */

public class TaskPool {

    static class TaskIndices {

        final int start;
        final int stop;
        final int task;

        TaskIndices(int start, int stop, int task) {
            this.start = start;
            this.stop = stop;
            this.task = task;
        }

        public String toString() {
            return start + " " + stop;
        }
    }

    private ExecutorService pool = null;
    private volatile PersistentForkCoordinator persistentForkCoordinator = null;
    private volatile DynamicForkCoordinator dynamicForkCoordinator = null;
    private volatile BalancedDynamicForkCoordinator balancedDynamicForkCoordinator = null;
    final private List<TaskIndices> indices;
    final private int taskCount;
    final private int threadCount;

    public TaskPool(int taskCount, int threadCount) {
        this.indices = setupTasks(taskCount, Math.abs(threadCount));
        this.taskCount = taskCount;
        this.threadCount = threadCount;
    }

    public ExecutorService getPool() { return pool; }

    public List<TaskIndices> getIndices() { return indices; }

    public int getNumThreads() { return indices.size(); }

    public int getNumTaxon() { return taskCount; }

    private List<TaskIndices> setupTasks(int taskCount, int threadCount) {
        List<TaskIndices> tasks = new ArrayList<>(threadCount);

        int length = taskCount / threadCount;
        if (taskCount % threadCount != 0) ++length;

        int start = 0;

        for (int task = 0; task < threadCount && start < taskCount; ++task) {
            tasks.add(new TaskIndices(start, Math.min(start + length, taskCount), task));
            start += length;
        }

        return tasks;
    }

    private ExecutorService setupParallelServices(int threadCount) {

        final ExecutorService pool;

        if (threadCount > 1) {
            pool = Executors.newFixedThreadPool(threadCount);
        } else if (threadCount < 0) {
            pool = Executors.newCachedThreadPool();
        } else {
            pool = null;
        }

        return pool;
    }

    public interface TaskCallable {
        void execute(int task, int thread);
    }

    public interface RangeCallable<E> {
        E map(int start, int end, int thread);
    }

    private final class PersistentForkCoordinator {

        private final Thread[] workers;
        private final AtomicInteger completedWorkers;
        private volatile TaskCallable currentRunnable;
        private volatile Throwable failure;
        private volatile Thread waitingThread;
        private volatile int generation;

        private PersistentForkCoordinator() {
            this.workers = new Thread[indices.size()];
            this.completedWorkers = new AtomicInteger(0);
            this.generation = 0;
            for (final TaskIndices indexSet : indices) {
                final Thread worker = new Thread(
                        () -> runWorker(indexSet),
                        "TaskPool-worker-" + System.identityHashCode(TaskPool.this) + "-" + indexSet.task);
                worker.setDaemon(true);
                worker.start();
                workers[indexSet.task] = worker;
            }
        }

        private void runWorker(final TaskIndices indexSet) {
            int observedGeneration = 0;
            while (true) {
                while (generation == observedGeneration) {
                    LockSupport.park(this);
                }
                observedGeneration = generation;

                final TaskCallable runnable = currentRunnable;
                if (runnable == null) {
                    return;
                }

                try {
                    for (int task = indexSet.start; task < indexSet.stop; ++task) {
                        runnable.execute(task, indexSet.task);
                    }
                } catch (Throwable throwable) {
                    recordFailure(throwable);
                }

                if (completedWorkers.incrementAndGet() == workers.length) {
                    LockSupport.unpark(waitingThread);
                }
            }
        }

        private synchronized void execute(final TaskCallable runnable) {
            failure = null;
            completedWorkers.set(0);
            waitingThread = Thread.currentThread();
            currentRunnable = runnable;
            generation++;
            for (final Thread worker : workers) {
                LockSupport.unpark(worker);
            }

            while (completedWorkers.get() < workers.length) {
                LockSupport.park(this);
            }
            waitingThread = null;

            final Throwable throwable = failure;
            if (throwable != null) {
                if (throwable instanceof RuntimeException) {
                    throw (RuntimeException) throwable;
                }
                if (throwable instanceof Error) {
                    throw (Error) throwable;
                }
                throw new RuntimeException(throwable);
            }
        }

        private void recordFailure(final Throwable throwable) {
            if (failure == null) {
                synchronized (this) {
                    if (failure == null) {
                        failure = throwable;
                    }
                }
            }
        }
    }

    private PersistentForkCoordinator getPersistentForkCoordinator() {
        PersistentForkCoordinator coordinator = persistentForkCoordinator;
        if (coordinator == null) {
            synchronized (this) {
                coordinator = persistentForkCoordinator;
                if (coordinator == null) {
                    coordinator = new PersistentForkCoordinator();
                    persistentForkCoordinator = coordinator;
                }
            }
        }
        return coordinator;
    }

    private final class DynamicForkCoordinator {

        private final Thread[] workers;
        private final AtomicInteger completedWorkers;
        private final AtomicInteger nextTask;
        private volatile TaskCallable currentRunnable;
        private volatile Throwable failure;
        private volatile Thread waitingThread;
        private volatile int generation;
        private volatile int taskLimit;
        private volatile int chunkSize;

        private DynamicForkCoordinator() {
            this.workers = new Thread[indices.size()];
            this.completedWorkers = new AtomicInteger(0);
            this.nextTask = new AtomicInteger(0);
            this.generation = 0;
            for (int workerIndex = 0; workerIndex < workers.length; ++workerIndex) {
                final int thread = workerIndex;
                final Thread worker = new Thread(
                        () -> runWorker(thread),
                        "TaskPool-dynamic-worker-" + System.identityHashCode(TaskPool.this) + "-" + thread);
                worker.setDaemon(true);
                worker.start();
                workers[workerIndex] = worker;
            }
        }

        private void runWorker(final int thread) {
            int observedGeneration = 0;
            while (true) {
                while (generation == observedGeneration) {
                    LockSupport.park(this);
                }
                observedGeneration = generation;

                final TaskCallable runnable = currentRunnable;
                if (runnable == null) {
                    return;
                }

                try {
                    while (true) {
                        final int start = nextTask.getAndAdd(chunkSize);
                        if (start >= taskLimit) {
                            break;
                        }

                        final int stop = Math.min(start + chunkSize, taskLimit);
                        for (int task = start; task < stop; ++task) {
                            runnable.execute(task, thread);
                        }
                    }
                } catch (Throwable throwable) {
                    recordFailure(throwable);
                }

                if (completedWorkers.incrementAndGet() == workers.length) {
                    LockSupport.unpark(waitingThread);
                }
            }
        }

        private synchronized void execute(final int taskLimit,
                                          final int chunkSize,
                                          final TaskCallable runnable) {
            failure = null;
            completedWorkers.set(0);
            nextTask.set(0);
            waitingThread = Thread.currentThread();
            currentRunnable = runnable;
            this.taskLimit = taskLimit;
            this.chunkSize = Math.max(1, chunkSize);
            generation++;
            for (final Thread worker : workers) {
                LockSupport.unpark(worker);
            }

            while (completedWorkers.get() < workers.length) {
                LockSupport.park(this);
            }
            waitingThread = null;

            final Throwable throwable = failure;
            if (throwable != null) {
                if (throwable instanceof RuntimeException) {
                    throw (RuntimeException) throwable;
                }
                if (throwable instanceof Error) {
                    throw (Error) throwable;
                }
                throw new RuntimeException(throwable);
            }
        }

        private void recordFailure(final Throwable throwable) {
            if (failure == null) {
                synchronized (this) {
                    if (failure == null) {
                        failure = throwable;
                    }
                }
            }
        }
    }

    private DynamicForkCoordinator getDynamicForkCoordinator() {
        DynamicForkCoordinator coordinator = dynamicForkCoordinator;
        if (coordinator == null) {
            synchronized (this) {
                coordinator = dynamicForkCoordinator;
                if (coordinator == null) {
                    coordinator = new DynamicForkCoordinator();
                    dynamicForkCoordinator = coordinator;
                }
            }
        }
        return coordinator;
    }

    private final class BalancedDynamicForkCoordinator {

        private final Thread[] workers;
        private final AtomicInteger completedWorkers;
        private final AtomicInteger nextTask;
        private volatile TaskCallable currentRunnable;
        private volatile Throwable failure;
        private volatile Thread waitingThread;
        private volatile int generation;
        private volatile int taskLimit;
        private volatile int chunkSize;
        private volatile int initialChunkOffset;

        private BalancedDynamicForkCoordinator() {
            this.workers = new Thread[indices.size()];
            this.completedWorkers = new AtomicInteger(0);
            this.nextTask = new AtomicInteger(0);
            this.generation = 0;
            this.initialChunkOffset = 0;
            for (int workerIndex = 0; workerIndex < workers.length; ++workerIndex) {
                final int thread = workerIndex;
                final Thread worker = new Thread(
                        () -> runWorker(thread),
                        "TaskPool-balanced-worker-" + System.identityHashCode(TaskPool.this) + "-" + thread);
                worker.setDaemon(true);
                worker.start();
                workers[workerIndex] = worker;
            }
        }

        private void runWorker(final int thread) {
            int observedGeneration = 0;
            while (true) {
                while (generation == observedGeneration) {
                    LockSupport.park(this);
                }
                observedGeneration = generation;

                final TaskCallable runnable = currentRunnable;
                if (runnable == null) {
                    return;
                }

                try {
                    runParticipant(thread, runnable);
                } catch (Throwable throwable) {
                    recordFailure(throwable);
                }

                if (completedWorkers.incrementAndGet() == workers.length) {
                    LockSupport.unpark(waitingThread);
                }
            }
        }

        private void runParticipant(final int thread, final TaskCallable runnable) {
            final int effectiveChunkSize = chunkSize;
            final int participantCount = workers.length + 1;
            final int initialChunkIndex = (thread + initialChunkOffset) % participantCount;
            final int initialStart = initialChunkIndex * effectiveChunkSize;

            if (initialStart < taskLimit) {
                final int initialStop = Math.min(initialStart + effectiveChunkSize, taskLimit);
                for (int task = initialStart; task < initialStop; ++task) {
                    runnable.execute(task, thread);
                }
            }

            while (true) {
                final int start = nextTask.getAndAdd(effectiveChunkSize);
                if (start >= taskLimit) {
                    break;
                }

                final int stop = Math.min(start + effectiveChunkSize, taskLimit);
                for (int task = start; task < stop; ++task) {
                    runnable.execute(task, thread);
                }
            }
        }

        private synchronized void execute(final int taskLimit,
                                          final int chunkSize,
                                          final TaskCallable runnable) {
            failure = null;
            completedWorkers.set(0);
            waitingThread = Thread.currentThread();
            currentRunnable = runnable;
            this.taskLimit = taskLimit;
            this.chunkSize = Math.max(1, chunkSize);
            final int participantCount = workers.length + 1;
            this.initialChunkOffset = generation % participantCount;
            nextTask.set(participantCount * this.chunkSize);
            generation++;
            for (final Thread worker : workers) {
                LockSupport.unpark(worker);
            }

            try {
                runParticipant(workers.length, runnable);
            } catch (Throwable throwable) {
                recordFailure(throwable);
            }

            while (completedWorkers.get() < workers.length) {
                LockSupport.park(this);
            }
            waitingThread = null;

            final Throwable throwable = failure;
            if (throwable != null) {
                if (throwable instanceof RuntimeException) {
                    throw (RuntimeException) throwable;
                }
                if (throwable instanceof Error) {
                    throw (Error) throwable;
                }
                throw new RuntimeException(throwable);
            }
        }

        private void recordFailure(final Throwable throwable) {
            if (failure == null) {
                synchronized (this) {
                    if (failure == null) {
                        failure = throwable;
                    }
                }
            }
        }
    }

    private BalancedDynamicForkCoordinator getBalancedDynamicForkCoordinator() {
        BalancedDynamicForkCoordinator coordinator = balancedDynamicForkCoordinator;
        if (coordinator == null) {
            synchronized (this) {
                coordinator = balancedDynamicForkCoordinator;
                if (coordinator == null) {
                    coordinator = new BalancedDynamicForkCoordinator();
                    balancedDynamicForkCoordinator = coordinator;
                }
            }
        }
        return coordinator;
    }

    public <E> E mapReduce(final RangeCallable<E> map, final BinaryOperator<E> reduce) {

        E result = null;

        if (indices.size() == 1) {

            final TaskIndices index = indices.get(0);
            result = map.map(index.start, index.stop, 0);

        } else {

            if (pool == null) {
                pool = setupParallelServices(threadCount);
            }

            List<Callable<E>> calls = new ArrayList<>();

            for (final TaskIndices indexSet : indices) {
                calls.add(() -> map.map(indexSet.start, indexSet.stop, indexSet.task));
            }

            try {

                List<Future<E>> futures = pool.invokeAll(calls);

                result = futures.get(0).get();
                for (int i = 1; i < futures.size(); ++i) {
                    result = reduce.apply(result, futures.get(i).get());
                }

            } catch (InterruptedException | ExecutionException exception) {
                exception.printStackTrace();
            }

        }

        return result;
    }

    public void fork(final TaskCallable runnable) {
        if (indices.size() == 1) {

            final TaskIndices index = indices.get(0);
            for (int task = index.start; task < index.stop; ++task) {
                runnable.execute(task, 0);
            }

        } else {
            if (threadCount > 1) {
                getPersistentForkCoordinator().execute(runnable);
                return;
            }

            if (pool == null) {
                pool = setupParallelServices(threadCount);
            }

            List<Callable<Object>> calls = new ArrayList<>();

            for (final TaskIndices indexSet : indices) {

                calls.add(Executors.callable(() -> {
                            for (int task = indexSet.start; task < indexSet.stop; ++task) {
                                runnable.execute(task, indexSet.task);
                            }
                        }
                ));
            }

            try {
                pool.invokeAll(calls);
            } catch (InterruptedException exception) {
                exception.printStackTrace();
            }
        }
    }

    public void forkDynamic(final int taskLimit,
                            final int chunkSize,
                            final TaskCallable runnable) {
        if (taskLimit <= 0) {
            return;
        }

        if (indices.size() == 1) {
            for (int task = 0; task < taskLimit; ++task) {
                runnable.execute(task, 0);
            }
        } else if (threadCount > 1) {
            getDynamicForkCoordinator().execute(taskLimit, chunkSize, runnable);
        } else {
            if (pool == null) {
                pool = setupParallelServices(threadCount);
            }

            final AtomicInteger nextTask = new AtomicInteger(0);
            final int effectiveChunkSize = Math.max(1, chunkSize);
            List<Callable<Object>> calls = new ArrayList<>();

            for (final TaskIndices indexSet : indices) {
                calls.add(Executors.callable(() -> {
                            while (true) {
                                final int start = nextTask.getAndAdd(effectiveChunkSize);
                                if (start >= taskLimit) {
                                    break;
                                }
                                final int stop = Math.min(start + effectiveChunkSize, taskLimit);
                                for (int task = start; task < stop; ++task) {
                                    runnable.execute(task, indexSet.task);
                                }
                            }
                        }
                ));
            }

            try {
                pool.invokeAll(calls);
            } catch (InterruptedException exception) {
                exception.printStackTrace();
            }
        }
    }

    public void forkDynamicBalanced(final int taskLimit,
                                    final int chunkSize,
                                    final TaskCallable runnable) {
        if (taskLimit <= 0) {
            return;
        }

        if (indices.size() == 1) {
            for (int task = 0; task < taskLimit; ++task) {
                runnable.execute(task, 0);
            }
        } else if (threadCount > 1) {
            getBalancedDynamicForkCoordinator().execute(taskLimit, chunkSize, runnable);
        } else {
            forkDynamic(taskLimit, chunkSize, runnable);
        }
    }
}
