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

import static dr.math.matrixAlgebra.missingData.InversionResult.Code.*;

/**
 * @author Marc A. Suchard
 */
public class InversionResult {

    public enum Code {
        FULLY_OBSERVED,
        NOT_OBSERVED,
        PARTIALLY_OBSERVED
    }

    InversionResult(Code code, int dim, double logDeterminant, boolean isLog) {
        this.code = code;
        this.dim = dim;
        this.logDeterminant = logDeterminant;
        this.isLog = isLog;
    }

    public InversionResult(Code code, int dim, double logDeterminant) {
        this(code, dim, logDeterminant, true);
    }

    final public Code getReturnCode() {
        return code;
    }

    final public int getEffectiveDimension() {
        assert (dim > -1) : "Should not try to get this effective dimension.";
        return dim;
    }

    final public double getDeterminant() {
        return isLog ? Math.exp(logDeterminant) : logDeterminant;
    }

    final public double getLogDeterminant() {
        return isLog ? logDeterminant : Math.log(logDeterminant);
    }

    public static InversionResult mult(InversionResult c1, InversionResult c2) {
        double logDet = c1.getLogDeterminant() + c2.getLogDeterminant();
        if (c1.getEffectiveDimension() == 0 || c2.getEffectiveDimension() == 0) {
            return new InversionResult(Code.NOT_OBSERVED, 0, logDet, true);
        }
        if (c1.getReturnCode() == Code.FULLY_OBSERVED) {
            return new InversionResult(c2.getReturnCode(), c2.getEffectiveDimension(), logDet, true);
        }
        if (c2.getReturnCode() == Code.FULLY_OBSERVED) {
            return new InversionResult(c1.getReturnCode(), c1.getEffectiveDimension(), logDet, true);
        }
        return new InversionResult(Code.PARTIALLY_OBSERVED, -1, logDet, true);
        // Effective dimension is unknown in this last case (<= min(dim1, dim2)), but should never be used.
    }

    public static Code getCode(int fullDim, int effectiveDim) {
        final InversionResult.Code code;
        if (effectiveDim == 0) {
            code = NOT_OBSERVED;
        } else if (effectiveDim == fullDim) {
            code = FULLY_OBSERVED;
        } else {
            code = PARTIALLY_OBSERVED;
        }
        return code;
    }

    public String toString() {
        return code + ":" + dim + ":" + logDeterminant;
    }

    final private Code code;
    final private int dim;
    final private double logDeterminant;
    final private boolean isLog;
}

