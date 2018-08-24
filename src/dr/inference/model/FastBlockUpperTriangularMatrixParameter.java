/*
 * BlockUpperTriangularMatrixParameter.java
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



public class FastBlockUpperTriangularMatrixParameter extends FastMatrixParameter {
    private int rowDim;
    private Bounds bounds = null;

//    public TransposedBlockUpperTriangularMatrixParameter transposeBlock(){
//        return TransposedBlockUpperTriangularMatrixParameter.recast(getVariableName(), this);
//    }

    public FastBlockUpperTriangularMatrixParameter(String name, int rows, int cols) {
        super(name, cols, rows, 1);


        boolean diagonalRestriction=false;
        if(diagonalRestriction) {
            for (int i = 0; i < getColumnDimension(); i++) {
                if (i < getRowDimension()) {
                    double[] uppers = new double[i + 1];
                    double[] lowers = new double[i + 1];

                    for (int j = 0; j < uppers.length; j++) {
                        uppers[j] = Double.POSITIVE_INFINITY;
                        lowers[j] = Double.NEGATIVE_INFINITY;

                    }
                    lowers[i] = 0;
                    getParameter(i).addBounds(new DefaultBounds(uppers, lowers));
                }
            }
        }

    }

//    @Override
//    public int getRowDimension() {
//        return rowDim;
//    }

//    public void setRowDimension(int rowDim){
//        this.rowDim=rowDim;
//    }

//    public double[][] getParameterAsMatrix(){
//        double[][] answer=new double[getRowDimension()][getColumnDimension()];
//        for(int i=0; i<getRowDimension(); i++){
//            for(int j=0; j<getColumnDimension(); j++){
//                if(i<=j){
////                    System.err.print(parameters[i].getSize());
////                    System.err.print("we get here\n");
//                    answer[i][j]=getParameter(j).getParameterValue(i);
//                }
//                else{
////                    System.err.print(i);
////                    System.err.print(" ");
////                    System.err.print(j);
////                    System.err.print("\n");
//                    answer[i][j]=0;
////                    System.err.print("getting here?\n");
//                }
//            }
//        }
//
//        return answer;
//    }

    public double getParameterValue(int row, int col) {
        if (!matrixCondition(row, col)) {
            return 0.0;
        } else {
            return super.getParameterValue(row,col);
        }
    }

    protected int getRow(int PID){
        return  PID%getRowDimension();
    }

    protected int getColumn(int PID){
        return PID/getRowDimension();
    }

    public void setParameterValueQuietly(int row, int col, double value){
        if(matrixCondition(row, col)){
            super.setParameterValueQuietly(row, col, value);
        }
    }

    public void setParameterValue(int row, int col,double value){
        if(matrixCondition(row, col)){
            super.setParameterValue(row, col, value);}
    }

    public void setParameterValue(int PID, double value){

        int row=getRow(PID);
        int col=getColumn(PID);
//        System.out.println(row+" "+col);
//        System.out.println(matrixCondition(row, col));


        if(matrixCondition(row, col)){
            setParameterValue(row, col, value);
        }
    }


    //test if violates matrix condition
    boolean matrixCondition(int row, int col){
        return row>=getColumnDimension()-row;
    }

    public double getParameterValue(int id){
        int row=getRow(id);
        int col=getColumn(id);

        if(matrixCondition(row, col)){
            return getParameterValue(row, col);
        }
        else
        {
            return 0;
        }
    }

    public void addBounds(Bounds<Double> boundary) {

        if (bounds == null) {
            bounds = new BUTMPBounds();
//            return;
        } //else {
        IntersectionBounds newBounds = new IntersectionBounds(getDimension());
        newBounds.addBounds(bounds);

//        }
        ((IntersectionBounds) bounds).addBounds(boundary);
    }

    public Bounds<Double> getBounds() {

        if (bounds == null) {
            bounds = new BUTMPBounds();
        }
        return bounds;
    }

    protected int getInnerDimension(int row, int col){
        return row;
    }


    private class BUTMPBounds implements Bounds<Double>{
        //TODO test!

        public Double getUpperLimit(int dim) {
            int row=getRow(dim);
            int col=getColumn(dim);

            if(matrixCondition(row, col)){

                return getParameter(col).getBounds().getUpperLimit(getInnerDimension(row, col)); }
            else
                return 0.0;
        }

        public Double getLowerLimit(int dim) {
            int row=getRow(dim);
            int col=getColumn(dim);

            if(matrixCondition(row, col)){
                return getParameter(col).getBounds().getLowerLimit(getInnerDimension(row, col));
//                    System.out.println(getParameters().get(dim-row).getBounds().getLowerLimit(getPindex().get(dim-row)));
//                return getParameters().get(dim-row).getBounds().getLowerLimit(getPindex().get(dim-row));
            }
            else
                return 0.0;
        }

        public int getBoundsDimension() {
//                int nBlanks = 0;
//                for (int i = 0; i <getColumnDimension() ; i++) {
//                    nBlanks+=1;
//                }
//
//                return getDimension()-nBlanks;
            return getDimension();
        }
    }

    public int getDimension(){
        return getRowDimension()*getColumnDimension();
    }


}
