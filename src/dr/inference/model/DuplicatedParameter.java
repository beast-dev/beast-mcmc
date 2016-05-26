/*
 * DuplicatedParameter.java
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
 * @author Marc Suchard
 */
public class DuplicatedParameter extends Parameter.Abstract implements VariableListener {

    public DuplicatedParameter(Parameter parameter) {
        this.parameter = parameter;
        parameter.addVariableListener(this);
        copies = 1;
        originalBounds = parameter.getBounds();
        bounds = originalBounds;
    }

    public void addDuplicationParameter(Parameter dupParameter) {
        this.dupParameter = dupParameter;
        dupParameter.addParameterListener(this);
        updateDuplication();
    }

    private void updateDuplication() {
        copies = (int) dupParameter.getParameterValue(0);
        final int originalLength = originalBounds.getBoundsDimension();
        double[] lowers = new double[getDimension()];
        double[] uppers = new double[getDimension()];
        for(int i=0; i<originalLength; i++) {
            lowers[i] = originalBounds.getLowerLimit(i);
            uppers[i] = originalBounds.getUpperLimit(i);
        }
        for(int i=1; i<copies; i++) {
            System.arraycopy(lowers,0,lowers,i*originalLength,originalLength);
            System.arraycopy(uppers,0,uppers,i*originalLength,originalLength);
        }
        bounds = new DefaultBounds(uppers,lowers);
    }

    public int getDimension() {
        return parameter.getDimension() * copies;
    }

    protected void storeValues() {
        parameter.storeParameterValues();
    }

    protected void restoreValues() {
        parameter.restoreParameterValues();
    }

    protected void acceptValues() {
        parameter.acceptParameterValues();
    }

    protected void adoptValues(Parameter source) {
        parameter.adoptParameterValues(source);
    }

    public double getParameterValue(int dim) {
        return parameter.getParameterValue(dim % parameter.getDimension());
    }

    public void setParameterValue(int dim, double value) {
        parameter.setParameterValue(dim % parameter.getDimension(), value);
        fireParameterChangedEvent(dim, Parameter.ChangeType.VALUE_CHANGED);
    }

    public void setParameterValueQuietly(int dim, double value) {
        parameter.setParameterValueQuietly(dim % parameter.getDimension(), value);
    }

    public void setParameterValueNotifyChangedAll(int dim, double value){
        parameter.setParameterValueNotifyChangedAll(dim % parameter.getDimension(), value);
    }

    public String getParameterName() {
        if (getId() == null)
            return "duplicated" + parameter.getParameterName();
        return getId();
    }

    public void addBounds(Bounds bounds) {
        throw new RuntimeException("Not yet implemented.");
    }

    public Bounds<Double> getBounds() {
        return bounds;
    }

    public void addDimension(int index, double value) {
        throw new RuntimeException("Not yet implemented.");
    }

    public double removeDimension(int index) {
        throw new RuntimeException("Not yet implemented.");
    }

    public void variableChangedEvent(Variable variable, int index, ChangeType type) {
        if (variable == dupParameter) {
            updateDuplication();
        }
        // Values have changed, so notify listeners
        for (int i = 0; i < copies; ++i) {
            // fire once for each duplication
            fireParameterChangedEvent(index + i * parameter.getDimension(), Parameter.ChangeType.VALUE_CHANGED);
        }
    }

    private final Parameter parameter;
    private Parameter dupParameter;
    private int copies;
    private Bounds<Double> bounds;
    private final Bounds<Double> originalBounds;

}
