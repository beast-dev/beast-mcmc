package dr.app.beauti.enumTypes;

/**
 * @author Alexei Drummond
 * @author Walter Xie
 */
public enum LocationSubstModelType {
    SYM_SUBST("Symmetric substitution model"),
        ASYM_SUBST("Asymmetric substitution model");

        LocationSubstModelType(String name) {
            this.name = name;
        }

        public String toString() {
            return name;
        }

        private final String name;
}
