/*
 * HashedMissingArray.java
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

package dr.evomodel.treedatalikelihood.continuous;

import dr.math.matrixAlgebra.Vector;

import java.util.Arrays;

public class HashedMissingArray {

    final private double[] array;

    public HashedMissingArray(final double[] array) {
        this.array = array;
    }

    public double[] getArray() {
        return array;
    }

    public double get(int index) {
        return array[index];
    }

    public int getLength() {
        return array.length;
    }

    @Override
    public int hashCode() {
        return Arrays.hashCode(array);
    }

    @Override
    public boolean equals(Object obj) {
        return obj instanceof HashedMissingArray && Arrays.equals(array,
                ((HashedMissingArray) obj).array);
    }

    public String toString() {
        return new Vector(array).toString();
    }
}
