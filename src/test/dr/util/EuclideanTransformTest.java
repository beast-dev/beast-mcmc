/*
 * EuclideanTransformTest.java
 *
 * Copyright (c) 2002-2019 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

import dr.util.EuclideanBallToRTransform;
import dr.util.Transform;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;
import test.dr.inference.trace.TraceCorrelationAssert;

import java.text.NumberFormat;
import java.util.Locale;

/**
 * @author Paul Bastide
 */

public class EuclideanTransformTest extends TraceCorrelationAssert {

    private int dim;
    private double[] unconstrained;
    private EuclideanBallToRTransform transform;

    private NumberFormat format = NumberFormat.getNumberInstance(Locale.ENGLISH);

    public EuclideanTransformTest(String name) {
        super(name);
    }

    public void setUp() throws Exception {
        super.setUp();

        format.setMaximumFractionDigits(5);

        dim = 6;

        transform = new EuclideanBallToRTransform(dim);

        unconstrained = new double[]{-1000, 0.03, 1000, 0.0, -0.00005, 1.0};
    }

    public void testTransformation() {
        System.out.println("\nTest LKJ transform.");

        double[] constrained = transform.inverse(unconstrained, 0, unconstrained.length);
        double[] inverseUnconstrained = transform.transform(constrained, 0, constrained.length);

        for (int k = 0; k < unconstrained.length; k++) {
            assertEquals("inverse transform k=" + k,
                    format.format(unconstrained[k]),
                    format.format(inverseUnconstrained[k]));
        }
    }

    public void testJacobian() {
        System.out.println("\nTest LKJ Jacobian.");

        // Matrix
        double[][] jacobianMat = transform.computeJacobianMatrixInverse(unconstrained);
        double jacobianDetBis = Math.log(CommonOps.det(new DenseMatrix64F(jacobianMat)));

        // Determinant
        double jacobianDet = (new Transform.InverseMultivariate(transform)).getLogJacobian(unconstrained, 0, unconstrained.length);

        System.out.println("Log Jacobiant Det direct=" + jacobianDet);
        System.out.println("Log Jacobiant Det matrix=" + jacobianDetBis);

        assertEquals("jacobian log det",
                format.format(jacobianDet),
                format.format(jacobianDetBis));
    }
}
