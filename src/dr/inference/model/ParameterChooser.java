/*
 * ParameterChooser.java
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
 * @author Joseph Heled
 *         Date: 31/08/2009
 */
public class ParameterChooser extends Variable.BaseNumerical<Double> implements Variable<Double>, ModelListener {
    private final ValuesPool pool;
    private final int which[];

    ParameterChooser(String name, ValuesPool pool, int which[]) {
        super(name);
        this.pool = pool;
        this.which = which;

        pool.addModelListener(this);
    }

    public int getSize() {
        return which.length;
    }

    public Double getValue(int index) {
        return pool.getValue(which[index]);
    }

    public void setValue(int index, Double value) {
        assert false;
    }

    public Double[] getValues() {
        final int size = getSize();
        Double[] copyOfValues = new Double[size];
        for (int i = 0; i < size; i++) {
            copyOfValues[i] = getValue(i) ;
        }
        return copyOfValues;
    }

    public void storeVariableValues() {
        //
    }

    public void restoreVariableValues() {
        //
    }

    public void acceptVariableValues() {
      //
    }

    private final Bounds<Double> bounds = new Bounds<Double>() {
        public Double getUpperLimit(int dimension) {
            return pool.getBounds().getUpperLimit(which[dimension]);
        }

        public Double getLowerLimit(int dimension) {
            return pool.getBounds().getLowerLimit(which[dimension]);
        }

        public int getBoundsDimension() {
            return getSize();
        }
    };

    public Bounds<Double> getBounds() {
        return bounds;
    }

    public void addBounds(Bounds<Double> bounds) {
        final Bounds<Double> b = getBounds();
        assert bounds.getBoundsDimension() == b.getBoundsDimension();

       for(int dim = 0; dim < bounds.getBoundsDimension(); ++dim) {
            if( bounds.getLowerLimit(dim) > b.getLowerLimit(dim) ||
                bounds.getUpperLimit(dim) < b.getUpperLimit(dim) ) {
                throw new RuntimeException("can't do that");
            }
       }     
    }

    public void modelChangedEvent(Model model, Object object, int index) {        
        for(int k = 0; k < which.length; ++k) {
            if( pool.hasChanged(which[k], object, index) ) {
                fireVariableChanged(k);
            }
        }
    }

    public void modelRestored(Model model) {
    }
}
