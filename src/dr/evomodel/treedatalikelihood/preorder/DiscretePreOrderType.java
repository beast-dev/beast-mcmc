package dr.evomodel.treedatalikelihood.preorder;

/**
 * @author Marc Suchard
 * @author Filippo Monti
 * @version $Id$
 */
public enum DiscretePreOrderType {

    BOTTOM(0, "bottom"),
    TOP(1, "top"),
    ROTATED(2, "rotated at top");

    DiscretePreOrderType(int type, String meaning) {
        this.type = type;
        this.meaning = meaning;
    }

    public int getType() {
        return type;
    }

    public String getMeaning() {
        return meaning;
    }

    private final int type;
    private final String meaning;
}
