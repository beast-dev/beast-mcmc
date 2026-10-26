/*
 * MaximizerUnweightedObjectiveTest.java
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

package test.dr.evomodel.treedatalikelihood;

import com.github.lbfgs4j.liblbfgs.Function;
import dr.evomodel.treedatalikelihood.discrete.MaximizerWrtParameter;
import dr.inference.distribution.DistributionLikelihood;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.Parameter;
import dr.math.distributions.GammaDistribution;
import dr.util.Transform;
import junit.framework.TestCase;

/**
 * By default the maximizer optimizes the unweighted log density (no transform Jacobian; what TreeLime uses), whose
 * optimum is the mode in parameter space. With includeJacobian it optimizes log L - logJacobian, HMC's transformed
 * density; its value used to carry the wrong Jacobian sign relative to its gradient (plan 3, B1).
 */
public class MaximizerUnweightedObjectiveTest extends TestCase {

    private static final double SHAPE = 3.0;
    private static final double SCALE = 0.5;

    public MaximizerUnweightedObjectiveTest(String name) {
        super(name);
    }

    public void testValueAndGradientAgreeUnderLogTransform() {
        for (boolean includeJacobian : new boolean[]{false, true}) {
            Parameter x = new Parameter.Default(new double[]{0.4, 2.2});
            Function function = maximizer(x, includeJacobian).getFunction();

            double[] y = new double[]{Math.log(0.4), Math.log(2.2)};
            double[] analytic = function.gradientAt(y.clone());
            final double h = 1E-6;
            for (int i = 0; i < y.length; i++) {
                double[] plus = y.clone();
                plus[i] += h;
                double[] minus = y.clone();
                minus[i] -= h;
                final double numeric = (function.valueAt(plus) - function.valueAt(minus)) / (2 * h);
                assertEquals("includeJacobian=" + includeJacobian + " dimension " + i, numeric, analytic[i],
                        1E-6 * Math.max(1.0, Math.abs(numeric)));
            }
        }
    }

    // the default, used for optimization: the mode in parameter space, (shape - 1) scale
    public void testDefaultOptimumIsTheModeInParameterSpace() {
        Parameter x = new Parameter.Default(new double[]{0.4, 2.2});
        maximizer(x, false).maximize();
        for (int i = 0; i < x.getDimension(); i++) {
            assertEquals((SHAPE - 1) * SCALE, x.getParameterValue(i), 1E-5);
        }
    }

    // with the Jacobian, the mode of the density of log x, shape * scale
    public void testJacobianOptimumIsTheModeOfTheTransformedDensity() {
        Parameter x = new Parameter.Default(new double[]{0.4, 2.2});
        maximizer(x, true).maximize();
        for (int i = 0; i < x.getDimension(); i++) {
            assertEquals(SHAPE * SCALE, x.getParameterValue(i), 1E-5);
        }
    }

    private static MaximizerWrtParameter maximizer(Parameter x, boolean includeJacobian) {
        GammaDistribution gamma = new GammaDistribution(SHAPE, SCALE);
        DistributionLikelihood likelihood = new DistributionLikelihood(gamma);
        likelihood.addData(x);
        GradientWrtParameterProvider gradient = new GradientWrtParameterProvider.ParameterWrapper(gamma, x, likelihood);
        Transform transform = new Transform.Array(Transform.LOG, x.getDimension(), null);
        return new MaximizerWrtParameter(likelihood, x, gradient, transform,
                new MaximizerWrtParameter.Settings(0, 1E-6, true, false, includeJacobian));
    }
}
