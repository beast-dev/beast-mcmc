/*
 * TransposedMatrixParameter.java
 *
 * Copyright (c) 2002-2012 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

import java.util.ArrayList;
import java.util.List;

/**
 * @author Max Tolkoff
 * @author Marc A. Suchard
 */
public class TransposedMatrixParameter extends MatrixParameter {

    public TransposedMatrixParameter(String name) {
        super(name + ".transpose");
    }

    public double getParameterValue(int row, int col) {
        return super.getParameterValue(col, row);
    }

    public int getRowDimension() {
        return getParameterCount();
    }

    public int getColumnDimension() {
        return getParameter(0).getDimension();
    }

    public Parameter getParameter(int index) {
        if (slices == null) {
            // construct vector_slices
            slices = new ArrayList<Parameter>();
            for (int i = 0; i < getColumnDimension(); ++i) {
                VectorSliceParameter thisSlice = new VectorSliceParameter("name", i);
                for (int j = 0; j < getRowDimension(); ++j) {
                    thisSlice.addParameter(super.getParameter(j));
                }
                slices.add(thisSlice);
            }
        }
        return slices.get(index);
    }

    public double getParameterValue(int dim) {
        throw new RuntimeException("Not yet implemented");
    }

    private List<Parameter> slices = null;

}
