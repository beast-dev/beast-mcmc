package dr.app.beauti.types;

/**
 * @author Alexei Drummond
 * @author Walter Xie
 */
public class MicroSatModelType {

    public enum RateProportionality {
        EQUAL_RATE("Equal Rate"),
        PROPORTIONAL_RATE("Proportional Rate"),
        ASYM_QUAD("Asymetric Quadratic Rate");

        RateProportionality(String displayName) {
            this.displayName = displayName;
        }

        public String toString() {
            return displayName;
        }

        private String displayName;
    }

    public enum MutationalBias {
        UNBIASED("Unbiased"),
        CONSTANT_BIAS("Constant Bias"),
        LINEAR_BIAS("Linear Bias");

        MutationalBias(String displayName) {
            this.displayName = displayName;
        }

        public String toString() {
            return displayName;
        }

        private String displayName;
    }

    public enum Phase {
        ONE_PHASE("One Phase"),
        TWO_PHASE("Two Phase"),
        TWO_PHASE_STAR("Two Phase *");

        Phase(String displayName) {
            this.displayName = displayName;
        }

        public String toString() {
            return displayName;
        }

        private String displayName;
    }

}
