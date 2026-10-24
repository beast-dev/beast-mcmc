/*
 * ShrinkageClockGradientTest.java
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

package test.dr.inference.distribution;

import dr.app.beast.BeastParser;
import dr.app.beast.BeastVersion;
import dr.inference.distribution.shrinkage.BayesianBridgeDistributionModel;
import dr.inference.distribution.shrinkage.JointBayesianBridgeDistributionModel;
import dr.inference.distribution.shrinkage.MarginalBayesianBridgeDistributionModel;
import dr.inference.model.Parameter;
import junit.framework.TestCase;

import java.io.StringReader;

/**
 * Shrinkage local clock (plan 3, E6): its gradient elements are switched off for now and must refuse with a
 * "not yet implemented" error; the Bayesian bridge's global-scale derivative behind them is checked against finite
 * differences so it is ready when they are switched on.
 */
public class ShrinkageClockGradientTest extends TestCase {

    private static final double[] INCREMENTS = {0.3, -0.7, 0.1, 0.9, -0.2, 0.45, -1.1, 0.05};

    public ShrinkageClockGradientTest(String name) {
        super(name);
    }

    public void testJointBridgeGlobalScaleDerivative() {
        Parameter localScale = new Parameter.Default(new double[]{10.0, 2.0, 5.0, 0.5, 3.0, 8.0, 1.5, 4.0});
        for (Parameter slabWidth : new Parameter[]{null, new Parameter.Default(2.0)}) {
            Parameter globalScale = new Parameter.Default(0.8);
            check(new JointBayesianBridgeDistributionModel(globalScale, localScale, new Parameter.Default(0.25),
                    slabWidth, INCREMENTS.length, false), globalScale, "joint slab=" + (slabWidth != null));
        }
    }

    public void testMarginalBridgeGlobalScaleDerivative() {
        for (double alpha : new double[]{1.0, 1.5, 0.5}) {
            Parameter globalScale = new Parameter.Default(0.7);
            check(new MarginalBayesianBridgeDistributionModel(globalScale, new Parameter.Default(alpha),
                    INCREMENTS.length, false), globalScale, "marginal alpha=" + alpha);
        }
    }

    public void testGlobalScaleGradientIsRefused() {
        String xml = "<beast>\n" +
                "  <taxa id=\"taxa\"><taxon id=\"a\"/><taxon id=\"b\"/><taxon id=\"c\"/></taxa>\n" +
                "  <newick id=\"startingTree\">((a:1,b:1):1,c:2);</newick>\n" +
                "  <treeModel id=\"treeModel\"><tree idref=\"startingTree\"/>\n" +
                "    <rootHeight><parameter id=\"treeModel.rootHeight\"/></rootHeight>\n" +
                "    <nodeHeights internalNodes=\"true\"><parameter id=\"treeModel.internalNodeHeights\"/></nodeHeights>\n" +
                "  </treeModel>\n" +
                "  <arbitraryBranchRates id=\"substBranchRates\" centerAtOne=\"false\"><treeModel idref=\"treeModel\"/>\n" +
                "    <rates><parameter id=\"substBranchRates.rates\" value=\"1.2 0.8 1.1 0.9\"/></rates>\n" +
                "  </arbitraryBranchRates>\n" +
                "  <bayesianBridgeDistribution id=\"bbDistribution\">\n" +
                "    <dimension><parameter idref=\"substBranchRates.rates\"/></dimension>\n" +
                "    <globalScale><parameter id=\"globalScale\" value=\"1.0\"/></globalScale>\n" +
                "    <exponent><parameter id=\"exponent\" value=\"1.0\"/></exponent>\n" +
                "  </bayesianBridgeDistribution>\n" +
                "  <autoCorrelatedRatesPrior id=\"ratesPrior\" scaling=\"none\" log=\"true\" operateOnIncrements=\"true\">\n" +
                "    <arbitraryBranchRates idref=\"substBranchRates\"/><bayesianBridgeDistribution idref=\"bbDistribution\"/>\n" +
                "  </autoCorrelatedRatesPrior>\n" +
                "  <bayesianBridgeGlobalScaleGradient id=\"gradient\"><autoCorrelatedRatesPrior idref=\"ratesPrior\"/>" +
                "</bayesianBridgeGlobalScaleGradient>\n" +
                "</beast>\n";
        BeastParser parser = new BeastParser(new String[]{}, null, false, false, false, new BeastVersion());
        try {
            parser.parse(new StringReader(xml), null);
            fail("the shrinkage-clock gradient must be refused for now");
        } catch (Exception e) {
            String message = e.toString();
            for (Throwable cause = e.getCause(); cause != null; cause = cause.getCause()) {
                message += " / " + cause;
            }
            assertTrue(message, message.contains("shrinkage local clock are not yet implemented"));
        }
    }

    private static void check(BayesianBridgeDistributionModel bridge, Parameter globalScale, String label) {
        final double analytic = bridge.getGradientLogDensityWrtGlobalScale(INCREMENTS);
        final double value = globalScale.getParameterValue(0);
        final double h = 1E-6;
        globalScale.setParameterValue(0, value + h);
        final double plus = bridge.logPdf(INCREMENTS);
        globalScale.setParameterValue(0, value - h);
        final double minus = bridge.logPdf(INCREMENTS);
        globalScale.setParameterValue(0, value);
        final double numeric = (plus - minus) / (2 * h);
        assertEquals(label, numeric, analytic, 1E-6 * Math.max(1.0, Math.abs(numeric)));
    }
}
