/*
 * SafeInvertibleMatrixParameter.java
 *
 * Copyright (c) 2002-2018 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.NormOps;

import static dr.math.matrixAlgebra.missingData.MissingOps.wrapSpherical;

public class EigenVectorsMatrix extends MatrixParameter {

    EigenVectorsMatrix(String name) {
        super(name + ".eigenVectors");
    }

    @Override
    public boolean check(){
        DenseMatrix64F M = wrapSpherical(getParameterValues(), 0, getColumnDimension());
        return NormOps.conditionP2(M) < 1.0E8;
    }
}
