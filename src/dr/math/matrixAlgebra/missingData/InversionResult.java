/*
 * InversionResult.java
 *
 * Copyright (c) 2002-2017 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.math.matrixAlgebra.missingData;

/**
 * @author Marc A. Suchard
 */
public class InversionResult {

    public enum Code {
        FULLY_OBSERVED,
        NOT_OBSERVED,
        PARTIALLY_OBSERVED
    }

    public InversionResult(Code code, int dim, double determinant) {
        this.code = code;
        this.dim = dim;
        this.determinant = determinant;
    }

    final public Code getReturnCode() {
        return code;
    }

    final public int getEffectiveDimension() {
        return dim;
    }

    final public double getDeterminant() {
        return determinant;
    }

    public String toString() {
        return code + ":" + dim + ":" + determinant;
    }

    final private Code code;
    final private int dim;
    final private double determinant;
}

