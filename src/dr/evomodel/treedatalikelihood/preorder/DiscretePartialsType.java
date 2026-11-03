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
            return BeaglePartialsType.TOP;
        }
    },
    BOTTOM_SPECTRAL(1 << 2, "bottom-spectral") {
        @Override
        public BeaglePartialsType getBeagleType() {
            return BeaglePartialsType.BOTTOM;
        }
    },
    TOP_SPECTRAL(1 << 3, "top-spectral") {
        @Override
        public BeaglePartialsType getBeagleType() {
            return BeaglePartialsType.TOP;
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

    public static DiscretePartialsType parse(String string) {
        for (DiscretePartialsType type : DiscretePartialsType.values()) {
            if (string.compareToIgnoreCase(type.getMeaning()) == 0) {
                return type;
            }
        }
        throw new IllegalArgumentException("Unknown discrete partials type");
    }

    private final int type;
    private final String meaning;
}
