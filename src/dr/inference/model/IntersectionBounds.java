/*
 * IntersectionBounds.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
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
 *
 */

package dr.inference.model;

import java.util.ArrayList;


public class IntersectionBounds implements Bounds<Double> {

    IntersectionBounds(int dimension) {
        this.dimension = dimension;
    }

    public void addBounds(Bounds<Double> boundary) {
        if (boundary.getBoundsDimension() != dimension) {
            throw new IllegalArgumentException("Incorrect dimension of bounds, expected " +
                    dimension + " but received " + boundary.getBoundsDimension());
        }
        if (bounds == null) {
            bounds = new ArrayList<Bounds<Double>>();
        }
        bounds.add(boundary);
    }

    /**
     * Gets the maximum lower limit of this parameter and all its slave parameters
     */
    public Double getLowerLimit(int index) {

        double lower = Double.NEGATIVE_INFINITY;
        if (bounds != null) {
            for (Bounds<Double> boundary : bounds) {
                if (boundary.getLowerLimit(index) > lower) {
                    lower = boundary.getLowerLimit(index);
                }
            }
        }
        return lower;
    }

    /**
     * Gets the minimum upper limit of this parameter and all its slave parameters
     */
    public Double getUpperLimit(int index) {

        double upper = Double.POSITIVE_INFINITY;
        if (bounds != null) {
            for (Bounds<Double> boundary : bounds) {
                if (boundary.getUpperLimit(index) < upper) {
                    upper = boundary.getUpperLimit(index);
                }
            }
        }
        return upper;
    }

    public int getBoundsDimension() {
        return dimension;
    }

    public String toString() {
        StringBuilder str = new StringBuilder("upper=[" + getUpperLimit(0));
        for (int i = 1; i < getBoundsDimension(); i++) {
            str.append(", ").append(getUpperLimit(i));
        }
        str.append("] lower=[").append(getLowerLimit(0));
        for (int i = 1; i < getBoundsDimension(); i++) {
            str.append(", ").append(getLowerLimit(i));
        }

        str.append("]");
        return str.toString();
    }

    private ArrayList<Bounds<Double>> bounds = null;
    private final int dimension;
}
