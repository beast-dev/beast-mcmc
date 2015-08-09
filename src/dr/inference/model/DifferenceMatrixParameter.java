/*
 * DifferenceMatrixParameter.java
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

import java.util.List;

import dr.inference.model.Bounds;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.inference.model.VariableListener;

public class DifferenceMatrixParameter extends MatrixParameter implements
        VariableListener{

    private Parameter parameter1;
    private Parameter parameter2;
    private Bounds bounds = null;

    public DifferenceMatrixParameter(Parameter parameter1, Parameter parameter2) {
        super(null);

        this.parameter1 = parameter1;
        this.parameter2 = parameter2;

        this.parameter1.addVariableListener(this);
        this.parameter2.addVariableListener(this);

    }// END: Constructor

    @Override
    public double getParameterValue(int dim) {
        double value=0;
        if(dim<Math.max(parameter1.getDimension(),parameter2.getDimension()))
        {
            if(parameter1.getDimension()==parameter2.getDimension())  {
                value = parameter1.getParameterValue(dim)
                        - parameter2.getParameterValue(dim);                   }
            else if(parameter1.getDimension()<parameter2.getDimension()){
                value=parameter1.getParameterValue(dim%parameter1.getDimension())-parameter2.getParameterValue(dim);

            }
            else
                value=parameter1.getParameterValue(dim)-parameter2.getParameterValue(dim% parameter2.getDimension());
        }
        else {
            System.out.println(dim);
            throw new RuntimeException("Index out of bounds.");
        }
        return value;
    }// END: getParameterValue

    public double getParameterValue(int row, int col){
        double val1;
        double val2;
        if(parameter1 instanceof MatrixParameter){
            val1=((MatrixParameter) parameter1).getParameterValue(row, col);
        }
        else{
            val1=getParameterValue(row);
        }
        if(parameter2 instanceof MatrixParameter){
            val2=((MatrixParameter) parameter2).getParameterValue(row, col);
        }
        else{
            val2=getParameterValue(row);
        }
        return val1-val2;
    }

    public int getParameterCount(){
        int pcount1=1;
        int pcount2=1;
        if(parameter1 instanceof MatrixParameter){
            pcount1=((MatrixParameter) parameter1).getParameterCount();
        }
        if(parameter2 instanceof MatrixParameter){
            pcount2=((MatrixParameter) parameter2).getParameterCount();
        }
        return Math.max(pcount1, pcount2);
    }

    @Override
    public Parameter getParameter(int index) {
        Parameter tempParam1;
        Parameter tempParam2;
        if(parameter1 instanceof CompoundParameter){
            tempParam1=((CompoundParameter) parameter1).getParameter(index);
        }
        else{
            tempParam1=parameter1;
        }
        if(parameter2 instanceof CompoundParameter){
            tempParam2=((CompoundParameter) parameter2).getParameter(index);
        }
        else{
            tempParam2=parameter2;
        }
        return new DifferenceMatrixParameter(tempParam1, tempParam2);
    }

    @Override
    public void setParameterValue(int dim, double value) {
        throw new RuntimeException("Not implemented");
    }

    @Override
    public void setParameterValueQuietly(int dim, double value) {
        throw new RuntimeException("Not implemented");
    }

    @Override
    public void setParameterValueNotifyChangedAll(int dim, double value) {
        throw new RuntimeException("Not implemented");
    }

//    @Override
//    public String getParameterName() {
//        if (getId() == null) {
//
//            StringBuilder sb = new StringBuilder("ratio");
//            sb.append(parameter1.getId()).append(".")
//                    .append(parameter2.getId());
//            setId(sb.toString());
//        }
//
//        return getId();
//    }// END: getParameterName


    @Override
    public int getDimension() {
        return Math.max(parameter1.getDimension(), parameter2.getDimension());
    }

    @Override
    public void addBounds(Bounds<Double> bounds) {
        this.bounds = bounds;
    }

    @Override
    public Bounds<Double> getBounds() {
        return bounds;
    }

    @Override
    public void addDimension(int index, double value) {
        throw new RuntimeException("Not yet implemented.");
    }

    @Override
    public double removeDimension(int index) {
        throw new RuntimeException("Not yet implemented.");
    }

    @Override
    public void variableChangedEvent(Variable variable, int index,
                                     dr.inference.model.Variable.ChangeType type) {
        fireParameterChangedEvent(index, type);
    }

    @Override
    protected void storeValues() {
        parameter1.storeParameterValues();
        parameter2.storeParameterValues();
    }

    @Override
    protected void restoreValues() {
        parameter1.restoreParameterValues();
        parameter2.restoreParameterValues();
    }

    @Override
    public int getRowDimension() {
        boolean param1 = parameter1 instanceof MatrixParameter;
        boolean param2 = parameter2 instanceof MatrixParameter;
        System.out.println("row");
        if (param1 & param2) {
            if (((MatrixParameter) parameter1).getRowDimension() == ((MatrixParameter) parameter2).getRowDimension()) {
                return ((MatrixParameter) parameter1).getRowDimension();
            } else {
                throw new RuntimeException("parameters not of the same length");
            }
        } else if (param1) {
            return ((MatrixParameter) parameter1).getRowDimension();
        } else if (param2) {
            return ((MatrixParameter) parameter2).getRowDimension();
        } else {
            return parameter1.getDimension();
        }
    }

    @Override
    public int getColumnDimension(){
        System.out.println("column");
        boolean param1=parameter1 instanceof MatrixParameter;
        boolean param2=parameter2 instanceof MatrixParameter;
        if(param1 & param2){
            if(((MatrixParameter)parameter1).getColumnDimension()==((MatrixParameter)parameter2).getColumnDimension()){
                return ((MatrixParameter)parameter1).getColumnDimension();
            }
            else{
                throw new RuntimeException("parameters not of the same length");
            }
        }
        else if(param1) {
            return ((MatrixParameter) parameter1).getColumnDimension();
        }
        else if(param2){
            return ((MatrixParameter)parameter2).getColumnDimension();
        }
        else{
            return 1;
        }
    }

//    @Override
//    protected void acceptValues() {
//        parameter1.acceptParameterValues();
//        parameter2.acceptParameterValues();
//    }
//
//    @Override
//    protected void adoptValues(Parameter source) {
//        throw new RuntimeException("Not implemented");
//    }

}// END: class
