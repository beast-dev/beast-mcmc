package dr.app.beauti.components;

import dr.app.beauti.components.continuous.ContinuousModelExtensionType;

/**
 * @author Gabriel Hassler
 * @author Marc Suchard
 */

public class GeneratorHelper {

    public static String extensionModelName(ContinuousModelExtensionType extensionType, String modelName) {
        String extendedName;
        switch (extensionType) {
            case RESIDUAL:
                extendedName = "residualModel";
                break;
            case LATENT_FACTORS:
                extendedName = "factorModel";
                break;
            case NONE:
                throw new IllegalArgumentException("Should not be called");
            default:
                throw new IllegalArgumentException("Unknown extension type");

        }

        return modelName + "." + extendedName;
    }

    public static String extensionPrecisionName(String modelName) {
        return modelName + ".extensionPrecision";
    }

    public static String extensionPrecisionPriorName(String modelName) {
        return modelName + ".extensionPrecisionPrior";
    }
}
