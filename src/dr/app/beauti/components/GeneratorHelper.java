package dr.app.beauti.components;

import dr.app.beauti.util.XMLWriter;
import dr.util.Attribute;

/**
 * @author Gabriel Hassler
 * @author Marc Suchard
 */

public class GeneratorHelper {

//    public static String extensionModelName(ContinuousModelExtensionType extensionType, String modelName) {
//        String extendedName;
//        switch (extensionType) {
//            case RESIDUAL:
//                extendedName = "residualModel";
//                break;
//            case LATENT_FACTORS:
//                extendedName = "factorModel";
//                break;
//            case NONE:
//                throw new IllegalArgumentException("Should not be called");
//            default:
//                throw new IllegalArgumentException("Unknown extension type");
//
//        }
//
//        return modelName + "." + extendedName;
//    }

//    public static String extensionPrecisionName(String modelName) {
//        return modelName + ".extensionPrecision";
//    }

//    public static String extensionPrecisionPriorName(String modelName) {
//        return modelName + ".extensionPrecisionPrior";
//    }


    public static void writeMatrixParameter(XMLWriter writer, String id, int p) {

        double[][] values = new double[p][p];
        for (int i = 0; i < p; i++) {
            values[i][i] = 1;
        }

        int rowDim = p;
        int colDim = p;


        writer.writeOpenTag("matrixParameter", new Attribute[]{
                new Attribute.Default<>("id", id)
        });

        for (int i = 0; i < rowDim; i++) {
            StringBuilder sb = new StringBuilder();

            for (int j = 0; j < colDim; j++) {
                if (j > 0) {
                    sb.append(" ");
                }

                sb.append(values[i][j]);
            }


            writer.writeTag("parameter",
                    new Attribute[]{
                            new Attribute.Default<>("value", sb.toString())
                    }, true);
        }
        writer.writeCloseTag("matrixParameter");
    }

}
