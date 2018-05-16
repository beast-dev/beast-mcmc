/*
 * LKJTransformTest.java
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

package test.dr.util;

import dr.math.matrixAlgebra.IllegalDimension;
import dr.math.matrixAlgebra.Matrix;
import dr.math.matrixAlgebra.SymmetricMatrix;
import dr.util.CorrelationToCholesky;
import dr.util.LKJCholeskyTransformConstrained;
import dr.util.LKJTransformConstrained;
import dr.util.Transform;
import test.dr.inference.trace.TraceCorrelationAssert;

import java.text.NumberFormat;
import java.util.Locale;

import static dr.math.matrixAlgebra.SymmetricMatrix.compoundCorrelationSymmetricMatrix;

/**
 * @author Paul Bastide
 */

public class LKJTransformTest extends TraceCorrelationAssert {

    private int dim;
    private double[] CPCs;
    private double[] CPCsLimit;
    private double[][] CPCsMatrix;
    private LKJTransformConstrained transform;
    private LKJCholeskyTransformConstrained transformChol;
    private CorrelationToCholesky transformCorrToChol;

    private NumberFormat format = NumberFormat.getNumberInstance(Locale.ENGLISH);

    public LKJTransformTest(String name) {
        super(name);
    }

    public void setUp() throws Exception {
        super.setUp();

        format.setMaximumFractionDigits(5);

        dim = 6;

        transform = new LKJTransformConstrained(dim);
        transformChol = new LKJCholeskyTransformConstrained(dim);
        transformCorrToChol = new CorrelationToCholesky(dim);

        CPCs = new double[]{0.12, -0.13, 0.14, -0.15, 0.16,
                -0.23, 0.24, -0.25, 0.26,
                0.34, -0.35, 0.36,
                -0.45, 0.46,
                0.56};

        CPCsMatrix = new double[dim][dim];
        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) {
                CPCsMatrix[i][j] = (2 * (j % 2) - 1) * ((i + 1) * 0.1 + (j + 1) * 0.01);
            }
        }

        CPCsLimit = new double[]{0.99, -0.99, 0.99, -0.99, 0.99,
                -0.99, 0.99, -0.99, 0.99,
                -0.99, 0.99, -0.99,
                0.99, -0.99,
                0.99};
    }

    public void testGetter() {
        System.out.println("\nTest getters.");

        int k = 0;
        for (int i = 0; i < dim - 1; i++) {
            for (int j = i + 1; j < dim; j++) {
                assertEquals("getter i=" + i + "j=" + j,
                        format.format(CPCs[k]),
                        format.format(CPCsMatrix[i][j]));
                k++;
            }
        }
    }

    public void testTransformation() {
        System.out.println("\nTest LKJ transform.");

        double[] transformedValue = transform.inverse(CPCs);
        double[] inverseTransformedValues = transform.transform(transformedValue);

        SymmetricMatrix R = compoundCorrelationSymmetricMatrix(transformedValue, dim);

        System.out.println("transformedValue=" + R);
        try {
            assertTrue("Positive Definite", R.isPD());
        } catch (IllegalDimension illegalDimension) {
            illegalDimension.printStackTrace();
        }

        System.out.println("iCPC=" + new Matrix(inverseTransformedValues, dim * (dim - 1) / 2, 1));

        assertEquals("size CPCs",
                format.format(CPCs.length),
                format.format(inverseTransformedValues.length));

        for (int k = 0; k < CPCs.length; k++) {
            assertEquals("inverse transform k=" + k,
                    format.format(CPCs[k]),
                    format.format(inverseTransformedValues[k]));
        }
    }

    public void testTransformationLimit() {
        System.out.println("\nTest LKJ transform on the border.");

        double[] transformedValue = transform.inverse(CPCsLimit);
        double[] inverseTransformedValues = transform.transform(transformedValue);

        SymmetricMatrix R = compoundCorrelationSymmetricMatrix(transformedValue, dim);
        ;
        System.out.println("transformedValue=" + R);
        try {
            assertTrue("Positive Definite", R.isPD());
        } catch (IllegalDimension illegalDimension) {
            illegalDimension.printStackTrace();
        }

        System.out.println("iCPC=" + new Matrix(inverseTransformedValues, dim * (dim - 1) / 2, 1));

        assertEquals("size CPCs",
                format.format(CPCsLimit.length),
                format.format(inverseTransformedValues.length));

        for (int k = 0; k < CPCsLimit.length; k++) {
            assertEquals("inverse transform k=" + k,
                    format.format(CPCsLimit[k]),
                    format.format(inverseTransformedValues[k]));
        }
    }

    public void testTransformationRecursion() {
        System.out.println("\nTest LKJ transform.");

        double[] transformedValue = transform.inverseRecursion(CPCs, 0, CPCs.length);
        double[] transformedValueChol = transform.inverse(CPCs);

        for (int k = 0; k < transformedValueChol.length; k++) {
            assertEquals("transform chol rec k=" + k,
                    format.format(transformedValueChol[k]),
                    format.format(transformedValue[k]));
        }

        double[] inverseTransformedValues = transform.transformRecursion(transformedValue, 0, CPCs.length);

        SymmetricMatrix R = compoundCorrelationSymmetricMatrix(transformedValue, dim);
        ;
        System.out.println("transformedValue=" + R);
        try {
            assertTrue("Positive Definite", R.isPD());
        } catch (IllegalDimension illegalDimension) {
            illegalDimension.printStackTrace();
        }

        System.out.println("iCPC=" + new Matrix(inverseTransformedValues, dim * (dim - 1) / 2, 1));

        assertEquals("size CPCs",
                format.format(CPCs.length),
                format.format(inverseTransformedValues.length));

        for (int k = 0; k < CPCs.length; k++) {
            assertEquals("inverse transform k=" + k,
                    format.format(CPCs[k]),
                    format.format(inverseTransformedValues[k]));
        }

    }

    public void testTransformationLimitRecursion() {
        System.out.println("\nTest LKJ transform on the border.");

        double[] transformedValue = transform.inverseRecursion(CPCsLimit, 0, CPCsLimit.length);
        double[] transformedValueChol = transform.inverse(CPCsLimit);

        for (int k = 0; k < transformedValueChol.length; k++) {
            assertEquals("transform chol rec k=" + k,
                    format.format(transformedValueChol[k]),
                    format.format(transformedValue[k]));
        }

        double[] inverseTransformedValues = transform.transformRecursion(transformedValue, 0, CPCsLimit.length);

        SymmetricMatrix R = compoundCorrelationSymmetricMatrix(transformedValue, dim);
        ;
        System.out.println("transformedValue=" + R);
        try {
            assertTrue("Positive Definite", R.isPD());
        } catch (IllegalDimension illegalDimension) {
            illegalDimension.printStackTrace();
        }

        System.out.println("iCPC=" + new Matrix(inverseTransformedValues, dim * (dim - 1) / 2, 1));

        assertEquals("size CPCs",
                format.format(CPCsLimit.length),
                format.format(inverseTransformedValues.length));

        for (int k = 0; k < CPCsLimit.length; k++) {
            assertEquals("inverse transform k=" + k,
                    format.format(CPCsLimit[k]),
                    format.format(inverseTransformedValues[k]));
        }
    }

    public void testJacobian() {
        System.out.println("\nTest LKJ Jacobian.");

        // Matrix
        double[][] jacobianMat = transform.computeJacobianMatrixInverse(CPCs);

        Matrix Jac = new Matrix(jacobianMat);
        System.out.println("Jacobian Matrix=" + Jac.transpose());

        assertEquals("size Jacobian Matrix",
                format.format(dim * (dim - 1) / 2),
                format.format(Jac.rows()));

        assertEquals("size Jacobian Matrix",
                format.format(dim * (dim - 1) / 2),
                format.format(Jac.columns()));

        // Determinant
        double jacobianDet = (new Transform.InverseMultivariable(transform)).getLogJacobian(CPCs, 0, CPCs.length);

        double jacobianDetBis = 0;
        for (int i = 0; i < jacobianMat[0].length; i++) {
            jacobianDetBis += Math.log(jacobianMat[i][i]);
        }

        System.out.println("Log Jacobiant Det direct=" + jacobianDet);
        System.out.println("Log Jacobiant Det matrix=" + jacobianDetBis);

        assertEquals("jacobian log det",
                format.format(jacobianDet),
                format.format(jacobianDetBis));
    }

    public void testJacobianCholesky() {
        System.out.println("\nTest LKJ Cholesky Jacobian.");

        // Matrix
        double[][] jacobianMat = transformChol.computeJacobianMatrixInverse(CPCs);

        Matrix Jac = new Matrix(jacobianMat);
        System.out.println("Jacobian Matrix=" + Jac.transpose());

        assertEquals("size Jacobian Matrix",
                format.format(dim * (dim - 1) / 2),
                format.format(Jac.rows()));

        assertEquals("size Jacobian Matrix",
                format.format(dim * (dim - 1) / 2),
                format.format(Jac.columns()));

        // Determinant
        double jacobianDet = (new Transform.InverseMultivariable(transformChol)).getLogJacobian(CPCs, 0, CPCs.length);

        double jacobianDetBis = 0;
        for (int i = 0; i < jacobianMat[0].length; i++) {
            jacobianDetBis += Math.log(jacobianMat[i][i]);
        }

        System.out.println("Log Jacobiant Det direct=" + jacobianDet);
        System.out.println("Log Jacobiant Det matrix=" + jacobianDetBis);

        assertEquals("jacobian log det",
                format.format(jacobianDet),
                format.format(jacobianDetBis));
    }

    public void testJacobianComposition() {
        System.out.println("\nTest LKJ Composition Cholesky Jacobian.");

        // Transforms
        double[] cholValues = transformChol.inverse(CPCs);
        double[] corrValues = transform.inverse(CPCs);
        double[] corrValuesBis = transformCorrToChol.inverse(cholValues);

        Transform.MultivariableTransform transformComposition = new Transform.ComposeMultivariable(transformChol, transformCorrToChol);
        double[] corrValuesTer = transformComposition.inverse(CPCs, 0, CPCs.length);


        for (int k = 0; k < CPCs.length; k++) {
            assertEquals("inverse transform k=" + k,
                    format.format(corrValues[k]),
                    format.format(corrValuesBis[k]));
        }

        for (int k = 0; k < CPCs.length; k++) {
            assertEquals("inverse transform k=" + k,
                    format.format(corrValues[k]),
                    format.format(corrValuesTer[k]));
        }

        // Determinant
        double jacobianDetCholToCPC = transformChol.getLogJacobian(cholValues, 0, CPCs.length);
        double jacobianDetCorrToChol = transformCorrToChol.getLogJacobian(corrValues, 0, CPCs.length);
        double jacobianDetCorrToCPC = transform.getLogJacobian(corrValues, 0, CPCs.length);
        double jacobianDetCorrToCPCComp = transformComposition.getLogJacobian(corrValues, 0, CPCs.length);

        System.out.println("Log Jacobiant Det Chol to CPC=" + jacobianDetCholToCPC);
        System.out.println("Log Jacobiant Det Corr to Chol=" + jacobianDetCorrToChol);
        System.out.println("Log Jacobiant Det Corr to CPC=" + jacobianDetCorrToCPC);
        System.out.println("Log Jacobiant Det Corr to CPC composition=" + jacobianDetCorrToCPCComp);

        assertEquals("jacobian log det",
                format.format(jacobianDetCorrToCPC),
                format.format(jacobianDetCholToCPC + jacobianDetCorrToChol));

        assertEquals("jacobian log det",
                format.format(jacobianDetCorrToCPC),
                format.format(jacobianDetCorrToCPCComp));
    }

}
