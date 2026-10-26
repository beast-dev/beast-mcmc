/*
 * GradientTypeTest.java
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

package test.dr.inference.hmc;

import dr.evomodel.treedatalikelihood.discrete.MaximizerWrtParameter;
import dr.inference.distribution.DistributionLikelihood;
import dr.inference.hmc.CompoundGradient;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.hmc.GradientWrtParameterProvider.GradientType;
import dr.inference.hmc.JointGradient;
import dr.inference.hmc.NumericalGradient;
import dr.inference.model.Parameter;
import dr.math.distributions.GammaDistribution;
import junit.framework.TestCase;

import java.util.Arrays;

/**
 * Gradient providers say whether they are analytic or numerical, and compound and joint gradients whether their
 * blocks are mixed; the maximizer reports it as "Gradient type".
 */
public class GradientTypeTest extends TestCase {

    public GradientTypeTest(String name) {
        super(name);
    }

    public void testBlocks() {
        Parameter x = new Parameter.Default("x", new double[]{0.4, 2.2});
        Parameter y = new Parameter.Default("y", new double[]{1.3});
        GradientWrtParameterProvider analyticX = analytic(x);
        GradientWrtParameterProvider numericalX = numerical(x);
        GradientWrtParameterProvider analyticY = analytic(y);
        GradientWrtParameterProvider numericalY = numerical(y);

        assertEquals(GradientType.ANALYTIC, analyticX.getGradientType());
        assertEquals(GradientType.NUMERICAL, numericalX.getGradientType());
        assertEquals(GradientType.NUMERICAL, new GradientWrtParameterProvider.Negative(numericalX).getGradientType());

        assertEquals(GradientType.ANALYTIC, new CompoundGradient(Arrays.asList(analyticX, analyticY)).getGradientType());
        assertEquals(GradientType.NUMERICAL, new CompoundGradient(Arrays.asList(numericalX, numericalY)).getGradientType());
        assertEquals(GradientType.MIXED, new CompoundGradient(Arrays.asList(analyticX, numericalY)).getGradientType());

        assertEquals(GradientType.ANALYTIC, new JointGradient(Arrays.asList(analyticX, analyticX)).getGradientType());
        assertEquals(GradientType.MIXED, new JointGradient(Arrays.asList(analyticX, numericalX)).getGradientType());
    }

    public void testMaximizerReport() {
        Parameter x = new Parameter.Default("x", new double[]{0.4, 2.2});
        Parameter y = new Parameter.Default("y", new double[]{1.3});
        CompoundGradient mixed = new CompoundGradient(Arrays.asList(analytic(x), numerical(y)));

        assertReports("numerical", numerical(x));
        assertReports("analytic", analytic(x));
        assertReports("mixed", mixed);

        // no gradient: the maximizer differences the density itself
        DistributionLikelihood likelihood = new DistributionLikelihood(new GammaDistribution(3.0, 0.5));
        likelihood.addData(x);
        MaximizerWrtParameter maximizer = new MaximizerWrtParameter(likelihood, x, null, null, settings());
        maximizer.maximize();
        assertTrue(maximizer.getReport().contains("Gradient type: numerical\n"));
    }

    private static void assertReports(String type, GradientWrtParameterProvider gradient) {
        MaximizerWrtParameter maximizer = new MaximizerWrtParameter(gradient.getLikelihood(), gradient.getParameter(),
                gradient, null, settings());
        maximizer.maximize();
        assertTrue(maximizer.getReport(), maximizer.getReport().contains("Gradient type: " + type + "\n"));
    }

    private static MaximizerWrtParameter.Settings settings() {
        return new MaximizerWrtParameter.Settings(3, 1E-6, true, false, false);
    }

    private static DistributionLikelihood gammaLikelihood(GammaDistribution gamma, Parameter x) {
        DistributionLikelihood likelihood = new DistributionLikelihood(gamma);
        likelihood.addData(x);
        return likelihood;
    }

    private static GradientWrtParameterProvider analytic(Parameter x) {
        GammaDistribution gamma = new GammaDistribution(3.0, 0.5);
        return new GradientWrtParameterProvider.ParameterWrapper(gamma, x, gammaLikelihood(gamma, x));
    }

    private static GradientWrtParameterProvider numerical(Parameter x) {
        return new NumericalGradient(gammaLikelihood(new GammaDistribution(3.0, 0.5), x), x);
    }
}
