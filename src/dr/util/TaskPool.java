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
import java.util.function.BinaryOperator;

/**
 * @author Marc A. Suchard
 * @author Andrew Holbrook
 */

public class TaskPool {

    class TaskIndices {

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
}
