/*
 * BayesianBridgeShrinkageOperatorTest.java
 *
 * Copyright (C) 2002-2009 Alexei Drummond and Andrew Rambaut
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
 * BEAST is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package test.dr.inference.operators;

import dr.inference.distribution.shrinkage.BayesianBridgeLikelihood;
import dr.inference.distribution.shrinkage.JointBayesianBridgeDistributionModel;
import dr.inference.model.Parameter;
import dr.inference.operators.shrinkage.BayesianBridgeShrinkageOperator;
import dr.math.MathUtils;
import dr.math.distributions.GammaDistribution;
import junit.framework.*;

public class BayesianBridgeShrinkageOperatorTest extends TestCase {

    public BayesianBridgeShrinkageOperatorTest(String name) {
        super(name);
    }

    public void setUp() throws Exception {
        super.setUp();
        MathUtils.setSeed(42);
    }

    public void testBayesianBridgeShrinkageOperator() {
        Parameter coefficients = new Parameter.Default("test.coefficents", 3, 0.0);
        Parameter exponent = new Parameter.Default("test.exponent", 1, 0.25);
        Parameter globalScale = new Parameter.Default("test.globalScale", 1, 0.25);
        Parameter localScale = new Parameter.Default("test.localScale", 3, 1.0);
        Parameter slabWidth = new Parameter.Default("test.slab", 1, 2.0);

        BayesianBridgeLikelihood bayesianBridge = new BayesianBridgeLikelihood(coefficients,
                new JointBayesianBridgeDistributionModel(globalScale, localScale, exponent, slabWidth,
                                                         coefficients.getDimension(), false));

        GammaDistribution globalScalePrior = new GammaDistribution(1.0, 1.0);

        BayesianBridgeShrinkageOperator bbso = new BayesianBridgeShrinkageOperator(bayesianBridge, globalScalePrior, null, 1.0);

        bbso.doOperation();
        double gs1 = globalScale.getParameterValue(0);
        double[] ls1 = localScale.getParameterValues();

        bbso.doOperation();
        double gs2 = globalScale.getParameterValue(0);
        double[] ls2 = localScale.getParameterValues();

        bbso.doOperation();
        double gs3 = globalScale.getParameterValue(0);
        double[] ls3 = localScale.getParameterValues();

//        System.err.println(gs1);
//        System.err.println(new dr.math.matrixAlgebra.Vector(ls1));
//        System.err.println(gs2);
//        System.err.println(new dr.math.matrixAlgebra.Vector(ls2));
//        System.err.println(gs3);
//        System.err.println(new dr.math.matrixAlgebra.Vector(ls3));

        assertEquals(gs1, 3.739465061574063E-5, 1e-8);
        assertEquals(ls1[0], 0.055002865424113125, 1e-8);
        assertEquals(ls1[1], 3.413792664116427E-7, 1e-8);
        assertEquals(ls1[2], 1.2209795305812487E-6, 1e-8);

        assertEquals(gs2, 7.459421891165917E-5, 1e-8);
        assertEquals(ls2[0], 2.042051812329108, 1e-8);
        assertEquals(ls2[1], 3.7658750918039714E-7, 1e-8);
        assertEquals(ls2[2], 1.40487675909568, 1e-8);

        assertEquals(gs3, 1.6334880920214307E-5, 1e-8);
        assertEquals(ls3[0], 0.005117931919904727, 1e-8);
        assertEquals(ls3[1], 0.7112807379133531, 1e-8);
        assertEquals(ls3[2], 69.80471269066086, 1e-8);

    }

    public static Test suite() {
        return new TestSuite(BayesianBridgeShrinkageOperatorTest.class);
    }
}