/*
 * TransposedBlockUpperTriangularMatrixParameter.java
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

import java.util.ArrayList;
import java.util.List;

/**
 * Created by max on 11/4/14.
 */
public class TransposedBlockUpperTriangularMatrixParameter extends BlockUpperTriangularMatrixParameter{
    public TransposedBlockUpperTriangularMatrixParameter(String name, Parameter[] params) {
        super(name, params, false);



        int colDim=params[params.length-1].getSize();
//        int rowDim=params.length;

//        for(int i=0; i<colDim; i++){
//            if(i<rowDim)
//            {params[i].setDimension(i+1);
//                this.addParameter(params[i]);}
//            else
//            {params[i].setDimension(rowDim);
//                this.addParameter(params[i]);
////                System.err.print(colDim-rowDim+i+1);
////                System.err.print("\n");
//            }
//        }
        this.colDim=colDim;
    }


    public static TransposedBlockUpperTriangularMatrixParameter recast(String name, CompoundParameter compoundParameter) {
        final int count = compoundParameter.getParameterCount();
        Parameter[] parameters = new Parameter[count];
        for (int i = 0; i < count; ++i) {
            parameters[i] = compoundParameter.getParameter(i);
        }
        return new TransposedBlockUpperTriangularMatrixParameter(name, parameters);
    }

//    public double getParameterValue(int row, int col){
//        if(col>row){
//            return 0;
//        }
//        else{
//            return getParameter(col).getParameterValue(row-col);
//        }
//    }

    protected int getRow(int PID){
        return  PID%getRowDimension();
    }

    protected int getColumn(int PID){
        return PID/getRowDimension();
    }

    @Override
    boolean matrixCondition(int row, int col) {
        return row>=col;
    }

    public void setParameterValue(int row, int col, double value){
        if(matrixCondition(row, col)){
            getParameter(col).setParameterValueQuietly(row - col, value);
            fireParameterChangedEvent(col*getRowDimension()+row, ChangeType.VALUE_CHANGED);
        }
    }

    public Parameter getParameter(int index) {
        if (slices == null) {
            // construct vector_slices
            slices = new ArrayList<Parameter>();
            for (int i = 0; i < getColumnDimension(); ++i) {
                VectorSliceParameter thisSlice = new VectorSliceParameter(getParameterName() + "." + i, i);
                for (int j = i; j < getRowDimension(); ++j) {
                    thisSlice.addParameter(super.getParameter(j));
                }
                slices.add(thisSlice);
            }
        }
        return slices.get(index);
    }

    protected int getInnerDimension(int row, int col){
        return row-col;
    }

    public int getRowDimension(){
        return getParameterCount();
    }

    public int getColumnDimension(){
        return colDim;
    }

    int colDim;
    private List<Parameter> slices = null;
}
