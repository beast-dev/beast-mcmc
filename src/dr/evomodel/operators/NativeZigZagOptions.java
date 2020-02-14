package dr.evomodel.operators;

/**
 * @author Marc A. Suchard
 */
public class NativeZigZagOptions {

    long flags;
    long seed;
    int info;

    public NativeZigZagOptions(NativeZigZag.Flag flags, long seed, int info) {
        this.flags = flags.getMask();
        this.seed = seed;
        this.info = info;
    }
}
