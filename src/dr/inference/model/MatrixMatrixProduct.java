/*
 * MatrixMatrixProduct.java
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

/*
@author Max Tolkoff
*/

//Designed to return a data matrix post computation if asked. Designed for latent liabilities
public class MatrixMatrixProduct extends MatrixParameter implements VariableListener {
    MatrixParameter left;
    MatrixParameter right;
    MatrixParameter inPlace;
    private final int leftDim;
    private final int rightDim;
    private final int midDim;
    Parameter columnMask;

    boolean[][] oldStoredValues;

    double[][] storedValues;
    boolean[][] areValuesStored;
    private Bounds bounds=null;

    public MatrixMatrixProduct(MatrixParameter[] params, Parameter columnMask) {
        super(null, params);



        this.columnMask=columnMask;

        this.left=params[0];
        this.right=params[1];
        if(params.length==3){
            inPlace=params[2];
            inPlace.addVariableListener(this);
        }
        storedValues=new double[left.getRowDimension()][right.getColumnDimension()];
        areValuesStored=new boolean[left.getRowDimension()][right.getColumnDimension()];
        oldStoredValues=new boolean[left.getRowDimension()][right.getColumnDimension()];

        for (int i = 0; i <left.getRowDimension() ; i++) {
            for (int j = 0; j <right.getColumnDimension() ; j++) {
                areValuesStored[i][j]=false;
            }

        }
        leftDim=left.getRowDimension();
        midDim=left.getColumnDimension();
        rightDim=right.getColumnDimension();

        left.addVariableListener(this);
        right.addVariableListener(this);
        inPlace.addVariableListener(this);

    }


    public void variableChangedEvent(Variable variable, int index, ChangeType type) {
        if(variable==right) {
//            System.out.println("RightChanged");
//            System.out.println(index/getRowDimension());
//            System.out.println(index);
            for (int i = 0; i <getRowDimension() ; i++) {
                areValuesStored[i][index/right.getRowDimension()]=false;
            }  }
        if(variable==left) {
//            System.out.println("LeftChanged");
//            System.out.println(index%left.getRowDimension());
//            System.out.println(index);
            for (int i = 0; i <getColumnDimension(); i++) {
                areValuesStored[index%left.getRowDimension()][i]=false;
            }  }
        fireParameterChangedEvent(index, type);
    }

    @Override
    public int getDimension(){
        return leftDim*rightDim;

    }

    public void addBounds(Bounds bounds) {
        this.bounds = bounds;
    }

    protected void storeValues() {
        System.arraycopy(areValuesStored, 0, oldStoredValues, 0, areValuesStored.length);
        left.storeParameterValues();
        right.storeParameterValues();
        inPlace.storeParameterValues();
    }

    protected void restoreValues() {
        left.restoreParameterValues();
        right.restoreVariableValues();
        inPlace.restoreParameterValues();
        areValuesStored=oldStoredValues;
    }

//    protected void acceptValues() {
//        left.acceptParameterValues();
//        right.acceptParameterValues();
//    }


    public double getParameterValue(int i, int j) {
        double sum = 0;
        if (columnMask.getParameterValue(j)!=0 && !areValuesStored[i][j]) {
            for (int k = 0; k < midDim; k++) {
                {
                    sum += left.getParameterValue(i, k) * right.getParameterValue(k, j);
                }

            }
        inPlace.setParameterValue(i,j, sum);
        areValuesStored[i][j]=true;
        }
        else{
            sum=inPlace.getParameterValue(i,j);
        }
        return sum;
    }

    @Override
    public double[][] getParameterAsMatrix() {
        return super.getParameterAsMatrix();
    }

    public Parameter getParameter(int PID) {
        for (int i = 0; i <leftDim ; i++) {
                    getParameterValue(i,PID);
                }
        return inPlace.getParameter(PID) ;
    }

    public int getRowDimension(){
        return leftDim;
    }

    public int getColumnDimension(){
        return rightDim;
    }

    @Override
    public double getParameterValue(int dim) {
        return getParameterValue(dim/getRowDimension(),dim%rightDim);
    }

    private void throwError(String functionName) throws RuntimeException {
        throw new RuntimeException("Object " + getId() + " is a deterministic function. Calling "
                + functionName + " is not allowed");
    }

    public void setParameterValue(int dim, double value) {
        throwError("setParameterValue()");
    }

    public void setParameterValueQuietly(int dim, double value) {
        throwError("setParameterValueQuietly()");
    }

    @Override

    public void setParameterValueNotifyChangedAll(int dim, double value) {
        throwError("setParameterValueNotifyChangedAll()");
    }


//    @Override
//    public String getParameterName() {
//        if (getId() == null) {
//            StringBuilder sb = new StringBuilder("product");
//            sb.append(".").append(left.getId());
//            sb.append(".").append(right.getId());
//            setId(sb.toString());
//        }
//        return getId();
//    }

    @Override
    public Bounds<Double> getBounds() {
        return bounds;
    }

    public void addDimension(int index, double value) {
        throw new RuntimeException("Not yet implemented.");
    }

    public double removeDimension(int index) {
        throw new RuntimeException("Not yet implemented.");
    }
};

