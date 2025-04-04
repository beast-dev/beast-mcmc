/*
 * GeneratorHelper.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 *  BEAST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 *
 */

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
