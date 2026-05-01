package dr.evomodel.treedatalikelihood.continuous.canonical;

/**
 * Strategy for choosing the chunk size used in dynamic task-pool parallelism.
 *
 * <p>A chunk size determines how many tasks each worker receives per scheduling
 * round. Larger chunks reduce scheduling overhead; smaller chunks improve load
 * balance when task durations vary.
 */
interface ChunkSizeStrategy {

    /**
     * Returns the chunk size to use when the task pool has {@code taskCount} tasks in total.
     *
     * @param taskCount total number of tasks to distribute (>= 1)
     * @return chunk size (>= 1)
     */
    int chunkSize(int taskCount);
}
