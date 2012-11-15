/*
 * ComplexColtEigenSystem.java
 *
 * Copyright (C) 2002-2012 Alexei Drummond, Andrew Rambaut & Marc A. Suchard
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

package dr.app.beagle.evomodel.substmodel;

import dr.math.matrixAlgebra.RobustEigenDecomposition;
import cern.colt.matrix.DoubleMatrix2D;

/**
 * @author Marc Suchard
 */
public class ComplexColtEigenSystem extends ColtEigenSystem {

    protected double[] getAllEigenValues(RobustEigenDecomposition decomposition) {
        double[] realEval = decomposition.getRealEigenvalues().toArray();
        double[] imagEval = decomposition.getImagEigenvalues().toArray();

        final int dim = realEval.length;
        double[] merge = new double[2*dim];
        System.arraycopy(realEval,0,merge,0,dim);
        System.arraycopy(imagEval,0,merge,dim,dim);
        return merge;
    }

      protected double[] getEmptyAllEigenValues(int dim) {
        return new double[2 * dim];
    }

    protected boolean validDecomposition(DoubleMatrix2D eigenV) {
        return true;
    }
}
