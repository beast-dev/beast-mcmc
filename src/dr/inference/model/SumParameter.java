/*
 * ProductParameter.java
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

import java.util.List;

/**
 * @author Marc Suchard
 */
public class SumParameter extends Parameter.Abstract implements VariableListener {

    public SumParameter(List<Parameter> parameterList) {
        this.parameterList = parameterList;
        dimension = parameterList.size() == 1 ? 1 : parameterList.get(0).getDimension();;
        for (Parameter p : parameterList) {
            p.addVariableListener(this);
        }
    }

    public int getDimension() {
        return dimension;
    }

    @Override
    public boolean isImmutable() {
        return true;
    }

    protected void storeValues() {
        for (Parameter p : parameterList) {
            p.storeParameterValues();
        }
    }

    protected void restoreValues() {
        for (Parameter p : parameterList) {
            p.restoreParameterValues();
        }
    }

    protected void acceptValues() {
        for (Parameter p : parameterList) {
            p.acceptParameterValues();
        }
    }

    protected void adoptValues(Parameter source) {
        throw new RuntimeException("Not implemented");
    }

    public double getParameterValue(int dim) {
        double value = 0;
        if (dimension == 1) {
            value = parameterList.get(0).getParameterValue(0);
            for (int i = 1; i < parameterList.get(0).getDimension(); i++) {
                value += parameterList.get(0).getParameterValue(i);
            }
        } else {
            value = parameterList.get(0).getParameterValue(dim);
            for (int i = 1; i < parameterList.size(); i++) {
                value += parameterList.get(i).getParameterValue(dim);
            }
        }
        return value;
    }

    public void setParameterValue(int dim, double value) {
        throw new RuntimeException("Not implemented");
    }

    public void setParameterValueQuietly(int dim, double value) {
        throw new RuntimeException("Not implemented");
    }

    public void setParameterValueNotifyChangedAll(int dim, double value){
        throw new RuntimeException("Not implemented");
    }

    public String getParameterName() {
        if (getId() == null) {
            StringBuilder sb = new StringBuilder("sum");
            for (Parameter p : parameterList) {
                sb.append(".").append(p.getId());
            }
            setId(sb.toString());
        }
        return getId();
    }

    public void addBounds(Bounds bounds) {
        this.bounds = bounds;
    }

    public Bounds<Double> getBounds() {
        if (bounds == null) {
            return parameterList.get(0).getBounds(); // TODO
        } else {
            return bounds;
        }
    }

    public void addDimension(int index, double value) {
        throw new RuntimeException("Not yet implemented.");
    }

    public double removeDimension(int index) {
        throw new RuntimeException("Not yet implemented.");
    }

    public void variableChangedEvent(Variable variable, int index, ChangeType type) {
        fireParameterChangedEvent(index,type);
    }

    private final List<Parameter> parameterList;
    private final int dimension;
    private Bounds bounds = null;
}
