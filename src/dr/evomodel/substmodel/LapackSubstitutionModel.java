/*
 * ComplexSubstitutionModel.java
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

import dr.evolution.datatype.DataType;
import dr.inference.model.BayesianStochasticSearchVariableSelection;
import dr.inference.model.Parameter;
import dr.math.matrixAlgebra.WrappedVector;
//import org.hipparchus.linear.MatrixUtils;
//import org.hipparchus.linear.RealMatrix;

//import static org.hipparchus.linear.MatrixUtils.createRealMatrix;

/**
 * @author Marc Suchard
 */

public class LapackSubstitutionModel extends ComplexSubstitutionModel {

    public LapackSubstitutionModel(String name, DataType dataType, FrequencyModel freqModel, Parameter parameter) {
        super(name, dataType, freqModel, parameter);
    }

    protected EigenSystem getDefaultEigenSystem(int stateCount) {
        return new LapackEigenSystem(stateCount);
    }

    private static final boolean CHECK_COLT = false;

    static double[] flatten(double[][] matrix) {

        int length = 0;
        for (double[] doubles : matrix) {
            length += doubles.length;
        }

        double[] result = new double[length];
        int offset = 0;
        for (double[] doubles : matrix) {
            System.arraycopy(doubles, 0, result, offset, doubles.length);
            offset += doubles.length;
        }

        return result;
    }

    public double getLogLikelihood() {

        if (CHECK_COLT) {

            EigenSystem ces = new ComplexColtEigenSystem(stateCount);
            EigenDecomposition ced = ces.decomposeMatrix(getQCopy());

            double[] cTest = new double[stateCount * stateCount];
            double[] lTest = new double[stateCount * stateCount];

            ces.computeExponential(ced, 1.0, cTest);

            if (cTest[0] == 0.0) {
                getTransitionProbabilities(1.0, lTest);

                double normalization = setupMatrix();
//                RealMatrix inf = createRealMatrix(getQCopy()).scalarMultiply(1 / normalization);
//                double[] gTest = flatten(MatrixUtils.matrixExponential(inf).getData());

                double[] q = flatten(getQCopy());
//                System.err.println("p   : " + new WrappedVector.Raw(gTest));
                System.err.println("p   : " + new WrappedVector.Raw(lTest));
                System.err.println("q   : " + new WrappedVector.Raw(q));
//                System.err.println("inf : " + new WrappedVector.Raw(flatten(inf.getData())));
                System.err.println("error");
            }
        }

        if (BayesianStochasticSearchVariableSelection.Utils.connectedAndWellConditioned(probability, this)) {
            return 0;
        }

        return Double.NEGATIVE_INFINITY;
    }
}
