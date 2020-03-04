package dr.evomodel.operators;

/**
 * @author Marc A. Suchard
 */
public class NativeZigZagOptions {

    private long flags;
    private long seed;
    private int info;

    public NativeZigZagOptions(long flags, long seed, int info) {
        this.flags = flags;
        this.seed = seed;
        this.info = info;
    }

    public long getFlags() { return flags; }

    public long getSeed() { return seed; }

    public int getInfo()  { return info; }
}
