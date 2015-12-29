/*
 * FiniteSetParameter.java
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
 * @author Marc A. Suchard
 */

public class FiniteSetParameter extends Parameter.Abstract implements VariableListener {

    public void variableChangedEvent(Variable variable, int index, ChangeType type) {
        //To change body of implemented methods use File | Settings | File Templates.
    }

    protected void storeValues() {
        indicator.storeParameterValues();
    }

    protected void restoreValues() {
        indicator.restoreParameterValues();
    }

    protected void acceptValues() {
        //To change body of implemented methods use File | Settings | File Templates.
    }

    protected void adoptValues(Parameter source) {
        //To change body of implemented methods use File | Settings | File Templates.
    }

    /**
     * @param dim the index of the parameter dimension of interest
     * @return the parameter's scalar value in the given dimension
     */
    public double getParameterValue(int dim) {
        return 0;  //To change body of implemented methods use File | Settings | File Templates.
    }

    /**
     * sets the scalar value in the given dimension of this parameter
     *
     * @param dim   the index of the dimension to set
     * @param value the value to set
     */
    public void setParameterValue(int dim, double value) {
        //To change body of implemented methods use File | Settings | File Templates.
    }

    /**
     * sets the scalar value in the given dimensin of this parameter to val, without firing any events
     *
     * @param dim   the index of the dimension to set
     * @param value the value to set
     */
    public void setParameterValueQuietly(int dim, double value) {
        //To change body of implemented methods use File | Settings | File Templates.
    }

    /**
     * sets the scalar value in the given dimensin of this parameter to val,
     * and notifies that values in all dimension have been changed
     *
     * @param dim   the index of the dimension to set
     * @param value the value to set
     */
    public void setParameterValueNotifyChangedAll(int dim, double value) {
        //To change body of implemented methods use File | Settings | File Templates.
    }

    /**
     * @return the name of this parameter
     */
    public String getParameterName() {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }

    /**
     * Adds new bounds to this parameter
     *
     * @param bounds to add
     */
    public void addBounds(Bounds<Double> bounds) {
        //To change body of implemented methods use File | Settings | File Templates.
    }

    /**
     * @return the intersection of all bounds added to this parameter
     */
    public Bounds<Double> getBounds() {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }

    /**
     * Adds an extra dimension at the given index
     *
     * @param index Index of the dimension to add
     * @param value value to save at end of new array
     */
    public void addDimension(int index, double value) {
        //To change body of implemented methods use File | Settings | File Templates.
    }

    /**
     * Removes the specified dimension from parameter
     *
     * @param index Index of dimension to lose
     * @return the value of the dimension removed
     */
    public double removeDimension(int index) {
        return 0;  //To change body of implemented methods use File | Settings | File Templates.
    }

    private Parameter indicator;

}
