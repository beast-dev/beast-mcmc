/*
 * EigenDecomposition.java
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

package dr.evomodel.substmodel;

//import org.hipparchus.linear.MatrixUtils;
//
//import static dr.evomodel.substmodel.LapackSubstitutionModel.flatten;

/**
 * @Author Marc A. Suchard
 */
public class LapackEigenDecomposition extends EigenDecomposition {

    final double[][] generator;

    private static final boolean CHECK_TRANSPOSE = false;

    public LapackEigenDecomposition(double[] evec, double[] ievc, double[] eval, double[][] generator) {
        super(evec, ievc, eval);
        this.generator = generator;
    }

    public LapackEigenDecomposition copy() {
        throw new RuntimeException("Not yet implemented");
    }

    public LapackEigenDecomposition transpose() {

        int dim = (int) Math.sqrt(Ievc.length);

        double[] test1;
//        double[] gTest1;

        if (CHECK_TRANSPOSE) {
            LapackEigenSystem les = new LapackEigenSystem(dim);
            test1 = new double[dim * dim];
            les.computeExponential(this, 1.0, test1);
//            gTest1 = flatten(MatrixUtils.matrixExponential(
//                    MatrixUtils.createRealMatrix(generator).scalarMultiply(1 / normalization)
//            ).getData());
            System.err.println("before");
        }

        EigenDecomposition tmp = super.transpose();

        LapackEigenDecomposition transposed = new LapackEigenDecomposition(tmp.Evec, tmp.Ievc, tmp.Eval, generator);

        double[] test2;
        double[] gTest2;

        if (CHECK_TRANSPOSE && Math.abs(test1[1] - test1[dim]) > 1E-3) { // prob[0][1] != prob[1][0]
            LapackEigenSystem les = new LapackEigenSystem(dim);
            test2 = new double[dim * dim];
            les.computeExponential(transposed, 1.0, test2);
//            gTest2 = gTest1.clone();
            gTest2 = test1.clone();
            transposeInPlace(gTest2, dim);
            System.err.println("after");
        }

        return transposed;
    }
}
