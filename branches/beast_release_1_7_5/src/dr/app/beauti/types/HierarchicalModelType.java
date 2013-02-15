package dr.app.beauti.types;

/**
 * @author Marc A. Suchard
 */
public enum HierarchicalModelType {
    
    NORMAL_HPM,
    LOGNORMAL_HPM;

    public String toString() {

        switch (this) {
            case NORMAL_HPM:
                return "Normal";
            case LOGNORMAL_HPM:
                return "Lognormal";            
            default:
                return "";
        }
    }
}
