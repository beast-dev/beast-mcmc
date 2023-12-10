package dr.evomodel.coalescent.basta;

import dr.evolution.tree.Tree;
import dr.math.matrixAlgebra.WrappedVector;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.*;

/**
 * @author Marc A. Suchard
 */
public class ParallelBastaLikelihoodDelegate extends GenericBastaLikelihoodDelegate {

    private static final int MIN_TASKS = 1;

    private final int threadCount;
    private final ExecutorService pool;
    private final List<Callable<Object>> allTasks;

    public ParallelBastaLikelihoodDelegate(String name,
                                           Tree tree,
                                           int stateCount,
                                           int threadCount) {
        super(name, tree, stateCount);

        if (threadCount > 1) {
            pool = Executors.newFixedThreadPool(threadCount);
        } else if (threadCount < -1) {
            pool = Executors.newCachedThreadPool();
        } else {
            throw new IllegalArgumentException("Illegal threadCount value");
        }

        this.threadCount = Math.abs(threadCount);

        this.allTasks = new ArrayList<>(threadCount);
        for (int i = 0; i < threadCount; ++i) {
            allTasks.add(null);
        }
    }

    protected void computeInnerBranchIntervalOperations(List<BranchIntervalOperation> branchIntervalOperations,
                                                      int start, int end) {

        int totalTasks = end - start;

        if (totalTasks <= MIN_TASKS) {
            super.computeInnerBranchIntervalOperations(branchIntervalOperations, start, end);
        } else {
            int numTasksPerThread = totalTasks / threadCount;
            if (totalTasks % threadCount != 0) {
                ++numTasksPerThread;
            }

            List<Callable<Object>> tasks = new ArrayList<>(threadCount);
            for (int i = 0, startTask = start; startTask < end; ++i, startTask += numTasksPerThread) {
                final int thisStart = startTask;
                final int thisEnd = Math.min(end, startTask + numTasksPerThread);

                tasks.add(() -> {
                    super.computeInnerBranchIntervalOperations(branchIntervalOperations, thisStart, thisEnd);
                    return null;
                });
            }
            try {
                pool.invokeAll(tasks);
            } catch (InterruptedException e) {
                throw new RuntimeException(e);
            }
        }
    }
}
