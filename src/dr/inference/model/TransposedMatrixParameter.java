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

    public TransposedMatrixParameter(String name, MatrixParameter param)
    {super(name + ".transpose");
        System.err.println("before");
        this.setRowDimension(param.getRowDimension());
        this.setColumnDimension(param.getColumnDimension());

        for(int i=0; i<param.getColumnDimension(); i++)
        {

            Parameter parameter=param.getParameter(i);
            System.err.println(param.getParameter(i).getDimension());
            this.addParameter(parameter);}
        System.err.println(param.getColumnDimension());
        System.err.println(param.getRowDimension());
        System.err.println(getRowDimension());
        System.err.println(getColumnDimension());
        System.err.println("middle");

        /*for(int i=1; i<getRowDimension(); i++){
            for(int j=0; j<getColumnDimension(); j++){
                this.addParameter(param.getParameter(j*getColumnDimension()+i));
            }}*/
     System.err.println("after");}

    /*public double getParameterValue(int row, int col) {
        return super.getParameterValue(col, row);
    }

    public int getRowDimension() {
        return getParameterCount();
    }

    public int getColumnDimension() {
        return getParameter(0).getDimension();
    }*/

    public void test()
    {
        System.err.println("I can read this!");
    }

    public Parameter getParameter(int index) {
        if (slices == null) {
            // construct vector_slices
            slices = new ArrayList<Parameter>();
            for (int i = 0; i < getColumnDimension(); ++i) {
                VectorSliceParameter thisSlice = new VectorSliceParameter("name", i);
                for (int j = 0; j < getRowDimension(); ++j) {
                    //does same thing each pass through i?
                    thisSlice.addParameter(super.getParameter(j));
                }
                slices.add(thisSlice);
            }
        }
        return slices.get(index);
    }

    /*public double getParameterValue(int row, int col) {
        return getParameter(row).getParameterValue(col);
    }

    public double[][] getParameterAsMatrix()
    {   final int J = getRowDimension();
        final int I = getColumnDimension();
        double[][] parameterAsMatrix = new double[I][J];
        for (int i = 0; i < I; i++) {
            for (int j = 0; j < J; j++)
                parameterAsMatrix[i][j] = getParameterValue(i, j);
        }
        return parameterAsMatrix;
    } */


    private List<Parameter> slices = null;

}
