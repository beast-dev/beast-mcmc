/*
 * MatrixParameterInterface.java
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

package dr.inference.model;

/**
 * Created by msuchard on 12/30/15.
 */
public interface MatrixParameterInterface extends Parameter {

    double getParameterValue(int row, int col);

    Parameter getParameter(int column); // Can return a proxy

    void setParameterValue(int row, int col, double value);

    void setParameterValueQuietly(int row, int col, double value);

    void setParameterValueNotifyChangedAll(int row, int col, double value);

    double[] getColumnValues(int col);

    double[][] getParameterAsMatrix();

    int getColumnDimension();

    int getRowDimension();

    double[] getParameterValues();

    int getUniqueParameterCount();

    Parameter getUniqueParameter(int index);

    void copyParameterValues(double[] destination, int offset);

//    void setAllParameterValuesQuietly(double[] values);

    void setAllParameterValuesQuietly(double[] values, int offset);

    String toSymmetricString();
}
