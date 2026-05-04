package dr.inference.timeseries.representation;

/**
 * Immutable snapshot of repeated-delta transition cache counters.
 */
public final class RepeatedDeltaCacheStatistics {

    public final int entryCount;
    public final long momentRequests;
    public final long momentBuilds;
    public final long canonicalRequests;
    public final long canonicalBuilds;

    RepeatedDeltaCacheStatistics(final int entryCount,
                                 final long momentRequests,
                                 final long momentBuilds,
                                 final long canonicalRequests,
                                 final long canonicalBuilds) {
        this.entryCount = entryCount;
        this.momentRequests = momentRequests;
        this.momentBuilds = momentBuilds;
        this.canonicalRequests = canonicalRequests;
        this.canonicalBuilds = canonicalBuilds;
    }

    public long momentHits() {
        return momentRequests - momentBuilds;
    }

    public long canonicalHits() {
        return canonicalRequests - canonicalBuilds;
    }
}
