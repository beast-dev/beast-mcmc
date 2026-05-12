package dr.evomodel.treedatalikelihood.preorder;

import beagle.BeaglePartialsType;

/**
 * @author Marc Suchard
 * @author Filippo Monti
 * @version $Id$
 */
public enum DiscretePartialsType {

    BOTTOM(1 << 0, "bottom") {
        @Override
        public BeaglePartialsType getBeagleType() {
            return BeaglePartialsType.BOTTOM;
        }
    },
    TOP(1 << 1, "top") {
        @Override
        public BeaglePartialsType getBeagleType() {
            return BeaglePartialsType.TOP.TOP;
        }
    };

    DiscretePartialsType(int type, String meaning) {
        this.type = type;
        this.meaning = meaning;
    }

    public int getType() {
        return type;
    }

    public String getMeaning() {
        return meaning;
    }

    abstract public BeaglePartialsType getBeagleType();

    private final int type;
    private final String meaning;
}
