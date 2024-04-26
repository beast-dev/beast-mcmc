package dr.app.beauti.util;

public class CharSetBlock {

    public CharSetBlock(int fromSite, int toSite, int every) {

        this.fromSite = fromSite;
        this.toSite = toSite;
        this.every = every;
    }

    public int getFromSite() {
        return fromSite;
    }

    public int getToSite() {
        return toSite;
    }

    public int getEvery() {
        return every;
    }

    private final int fromSite;
    private final int toSite;
    private final int every;
}
