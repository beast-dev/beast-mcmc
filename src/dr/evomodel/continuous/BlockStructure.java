package dr.evomodel.continuous;

public final class BlockStructure {

    private final int[] blockStarts;
    private final int[] blockSizes;

    public BlockStructure(int[] blockStarts, int[] blockSizes) {
        this.blockStarts = blockStarts.clone();
        this.blockSizes  = blockSizes.clone();
    }

    public int[] getBlockStarts() {
        return blockStarts.clone();
    }

    public int[] getBlockSizes() {
        return blockSizes.clone();
    }
}
