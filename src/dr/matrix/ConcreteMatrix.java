/*
 * ConcreteMatrix.java
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

package dr.matrix;

/**
 * @author Alexei Drummond
 */
class ConcreteMatrix extends MutableMatrix.AbstractMutableMatrix {

    public ConcreteMatrix(double[][] v) {

        this.values = new double[v.length][v[0].length];
        for (int i = 0; i < values.length; i++) {
            System.arraycopy(v[i], 0, values[i], 0, values[0].length);
        }
    }

    public final void setDimension(int rows, int columns) {

        if (values.length != rows || values[0].length != columns) {
            values = new double[rows][columns];
        }
    }

    public final void setElement(int rows, int column, double value) {
        values[rows][column] = value;
    }


    public final int getRowCount() {
        return values.length;
    }

    public final int getColumnCount() {
        return values[0].length;
    }

    public final double getElement(int i, int j) {
        return values[i][j];
    }

    public String toString() {

        StringBuffer buffer = new StringBuffer();
        for (int i = 0; i < values.length; i++) {
            for (int j = 0; j < values[0].length; j++) {
                buffer.append(values[i][j]).append("\t");
            }
            buffer.append("\n");
        }
        return buffer.toString();
    }

    double[][] values = null;
}
