package dr.evomodel.treedatalikelihood.continuous.canonical.scheduling;

/**
 * Chunk-size heuristic that scales with the trait dimension.
 *
 * <p>Higher-dimensional traits take more time per branch (O(d³) for matrix operations),
 * so fewer chunks per worker are needed to keep scheduling overhead low relative to
 * compute time. Conversely, small-d runs benefit from finer chunking to balance load.
 *
 * <p>Thresholds and limits were tuned empirically; they are exposed as named constants
 * so they can be documented, unit-tested, and adjusted in one place.
 */
public final class DimensionWeightedChunkSizeStrategy implements ChunkSizeStrategy {

    static final int HIGH_DIM_THRESHOLD   = 16;  // dim >= this → 1 chunk/worker, small max
    static final int MED_DIM_THRESHOLD    = 8;   // dim >= this → 2 chunks/worker, medium max
    static final int HIGH_DIM_CHUNKS_PER_WORKER = 1;
    static final int MED_DIM_CHUNKS_PER_WORKER  = 2;
    static final int LOW_DIM_CHUNKS_PER_WORKER  = 4;
    static final int HIGH_DIM_MAX_CHUNK   = 4;
    static final int MED_DIM_MAX_CHUNK    = 8;
    static final int LOW_DIM_MAX_CHUNK    = 32;

    private final int workerCount;
    private final int targetChunksPerWorker;
    private final int maxChunkSize;

    public DimensionWeightedChunkSizeStrategy(final int dimension, final int workerCount) {
        if (workerCount < 1) throw new IllegalArgumentException("workerCount must be >= 1");
        this.workerCount = workerCount;
        if (dimension >= HIGH_DIM_THRESHOLD) {
            this.targetChunksPerWorker = HIGH_DIM_CHUNKS_PER_WORKER;
            this.maxChunkSize          = HIGH_DIM_MAX_CHUNK;
        } else if (dimension >= MED_DIM_THRESHOLD) {
            this.targetChunksPerWorker = MED_DIM_CHUNKS_PER_WORKER;
            this.maxChunkSize          = MED_DIM_MAX_CHUNK;
        } else {
            this.targetChunksPerWorker = LOW_DIM_CHUNKS_PER_WORKER;
            this.maxChunkSize          = LOW_DIM_MAX_CHUNK;
        }
    }

    @Override
    public int chunkSize(final int taskCount) {
        final int suggested =
                (taskCount + workerCount * targetChunksPerWorker - 1)
                        / (workerCount * targetChunksPerWorker);
        return Math.max(1, Math.min(maxChunkSize, suggested));
    }
}
