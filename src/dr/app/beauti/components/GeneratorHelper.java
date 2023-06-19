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

    public static void writeParameter(XMLWriter writer, String id, double[] values) {


        if (values.length < 1) throw new IllegalArgumentException("cannot make parameter with < 1 dimensions");

        StringBuilder sb = new StringBuilder();
        sb.append(values[0]);
        for (int i = 1; i < values.length; i++) {
            sb.append(" ");
            sb.append(values[i]);
        }

        Attribute[] attributes;
        if (id == null) {
            attributes = new Attribute[]{
                    new Attribute.Default<>("value", sb.toString())
            };
        } else {
            attributes = new Attribute[]{
                    new Attribute.Default<String>("id", id),
                    new Attribute.Default<>("value", sb.toString())
            };
        }
        writer.writeTag("parameter", attributes, true);
    }

    public static void writeParameter(XMLWriter writer, String id, int p, double value) {
        double[] values = new double[p];
        for (int i = 0; i < p; i++) {
            values[i] = value;
        }
        writeParameter(writer, id, values);
    }

    public static void writeParameter(XMLWriter writer, String id, int p) {
        writeParameter(writer, id, p, 0);
    }


    public static void writeMatrixParameter(XMLWriter writer, String id, double[][] values) {

        int rowDim = values.length;

        writer.writeOpenTag("matrixParameter", new Attribute[]{
                new Attribute.Default<>("id", id)
        });

        for (int i = 0; i < rowDim; i++) {
            writeParameter(writer, null, values[i]);
        }
        writer.writeCloseTag("matrixParameter");
    }

    public static void writeMatrixParameter(XMLWriter writer, String id, int nRows, int nCols) {
        double[][] values = new double[nRows][nCols];
        writeMatrixParameter(writer, id, values);
    }

    public static void writeIdentityMatrixParameter(XMLWriter writer, String id, int p) {

        double[][] values = new double[p][p];
        for (int i = 0; i < p; i++) {
            values[i][i] = 1;
        }
        writeMatrixParameter(writer, id, values);
    }

    public static void writeMatrixInverse(XMLWriter writer, String id, String matrixId) {
        writer.writeOpenTag("matrixInverse",
                new Attribute[]{
                        new Attribute.Default<String>("id", id)
                });
        writer.writeIDref("matrixParameter", matrixId);
        writer.writeCloseTag("matrixInverse");
    }

}
