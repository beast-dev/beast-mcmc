/*
 * SubstitutionModelUtils.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
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
 */

package dr.oldevomodel.substmodel;

import dr.evolution.datatype.DataType;

import java.text.NumberFormat;

/**
 * @author Alexei Drummond
 */
public class SubstitutionModelUtils {

    public static String toString(double[] matrix, DataType dataType, boolean isTriangle, int indent) {

        int size = dataType.getStateCount();

        int expectedSize = size * size;
        if (isTriangle) expectedSize = (size * size) / 2 - size;

        assert matrix.length == expectedSize;

        double[][] amat = new double[size][size];

        int k = 0;
        for (int i = 0; i < size; i++) {
            for (int j = isTriangle ? i + 1 : 0; j < size; j++) {
                amat[i][j] = matrix[k];
                if (isTriangle) amat[j][i] = matrix[k];
                k += 1;
            }
        }
        return toString(amat, dataType, indent);
    }

    public static String frequenciesToString(double[] pi, DataType dataType, int indent) {

        int stateCount = pi.length;
        final int columnWidth = 7;

        StringBuilder builder = new StringBuilder();
        builder.append(header(dataType, indent, columnWidth));
        builder.append("\n");

        // write body
        NumberFormat formatter = NumberFormat.getNumberInstance();
        formatter.setMaximumFractionDigits(3);
        builder.append(spaces(indent));
        builder.append("  {");
        for (int j = 0; j < stateCount; j++) {

            builder.append(padded(formatter.format(pi[j]), columnWidth));
        }
        builder.append("}\n");
        return builder.toString();
    }


    public static String toString(double[][] matrix, DataType dataType, int indent) {

        int stateCount = matrix.length;
        final int columnWidth = 7;

        StringBuilder builder = new StringBuilder();
        builder.append(header(dataType, indent, columnWidth));
        builder.append("\n");

        // write matrix body

        NumberFormat formatter = NumberFormat.getNumberInstance();
        formatter.setMaximumFractionDigits(3);
        for (int i = 0; i < stateCount; i++) {
            builder.append(spaces(indent));
            builder.append(dataType.getChar(i)).append(" ");
            if (i == 0) {
                builder.append("/");
            } else if (i == stateCount - 1) {
                builder.append("\\");
            } else builder.append("|");
            for (int j = 0; j < stateCount; j++) {

                builder.append(padded(formatter.format(matrix[i][j]), columnWidth));
            }
            if (i == 0) {
                builder.append("\\\n");
            } else if (i == stateCount - 1) {
                builder.append("/\n");
            } else builder.append("|\n");
        }
        return builder.toString();
    }

    private static String spaces(int indent) {
        StringBuilder builder = new StringBuilder();
        for (int i = 0; i < indent; i++) {
            builder.append(' ');
        }
        return builder.toString();
    }

    private static String padded(String s, int width) {
        int extra = width - s.length();
        for (int i = 0; i < (extra / 2); i++) {
            s = " " + s;
        }
        extra = width - s.length();
        for (int i = 0; i < extra; i++) {
            s += " ";
        }
        return s;
    }

    private static String header(DataType dataType, int indent, int columnWidth) {
        // write header
        StringBuilder builder = new StringBuilder();
        builder.append(spaces(indent)).append("   ");
        for (int i = 0; i < dataType.getStateCount(); i++) {
            builder.append(padded(dataType.getChar(i) + "", columnWidth));
        }
        return builder.toString();
    }
}
