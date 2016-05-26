/*
 * Bounds.java
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

import java.io.Serializable;
import java.util.ArrayList;

/**
 * Represents a multi-dimensional 'regular' boundary (a hypervolume)
 *
 * @author Alexei Drummond
 * @version $Id: Bounds.java,v 1.2 2005/05/24 20:25:59 rambaut Exp $
 */
public interface Bounds<V> extends Serializable {

    /**
     * @return the upper limit of this hypervolume in the given dimension.
     */
    V getUpperLimit(int dimension);

    /**
     * @return the lower limit of this hypervolume in the given dimension.
     */
    V getLowerLimit(int dimension);

    /**
     * @return the dimensionality of this hypervolume.
     */
    int getBoundsDimension();

    public class Int implements Bounds<Integer> {

        int size = 1;
        int lower = java.lang.Integer.MIN_VALUE;
        int upper = java.lang.Integer.MAX_VALUE;

        public Int(int size, int lower, int upper) {
            this.size = size;
            this.lower = lower;
            this.upper = upper;
        }

        public Int(Variable<Integer> variable, int lower, int upper) {
            this.size = variable.getSize();
            this.lower = lower;
            this.upper = upper;
        }


        public Integer getUpperLimit(int dimension) {
            return upper;
        }

        public Integer getLowerLimit(int dimension) {
            return lower;
        }

        public int getBoundsDimension() {
            return size;
        }
    }

    /**
     * A staircase bound is used for model averaging and requires each index to have an upper bound equal to its
     * index. Thus there is always a value 0 in the first entry of the parameter, whereas there is a value of 0 or 1
     * in the second entry, {0,1,2} in the third entry et cetera. AJD
     *
     * But in the code, for simplicity, first entry of the parameter is not implemented,
     * so that integer parameter size = real size - 1, and 1st digital of paramter = 2nd position of proposed index,
     * namely {0, 1}.  Walter
     */
    public class Staircase implements Bounds<Integer> {

        int size = 0;

        public Staircase(int size) {
            this.size = size;
        }

        public Staircase(Variable<Integer> variable) {
            this.size = variable.getSize();
        }


        public Integer getUpperLimit(int dimension) {
            return dimension + 1; // integer index parameter size = real size - 1
        }

        public Integer getLowerLimit(int dimension) {
            return 0;
        }

        public int getBoundsDimension() {
            return size;
        }

        private ArrayList<Bounds<Integer>> bounds = null;

        public void addBounds(Bounds<Integer> boundary) {
            if (boundary.getBoundsDimension() != size) {
                throw new IllegalArgumentException("Incorrect dimension of bounds, expected " +
                        size + " but received " + boundary.getBoundsDimension());
            }
            if (bounds == null) {
                bounds = new ArrayList<Bounds<Integer>>();
            }
            bounds.add(boundary);
        }        
    }
}	
