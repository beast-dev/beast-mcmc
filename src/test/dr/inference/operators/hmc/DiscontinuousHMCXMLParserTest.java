/*
 * DiscontinuousHMCXMLParserTest.java
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

package test.dr.inference.operators.hmc;

import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.hmc.LikelihoodBasedDiscontinuousPotentialProvider;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.inference.model.ParameterParser;
import dr.inference.model.Variable;
import dr.inference.model.VariableListener;
import dr.inference.operators.AdaptationMode;
import dr.inference.operators.hmc.DiscontinuousHamiltonianMonteCarloOperator;
import dr.inference.operators.hmc.MixedDiscontinuousHamiltonianMonteCarloOperator;
import dr.inferencexml.hmc.LikelihoodBasedDiscontinuousPotentialProviderParser;
import dr.inferencexml.operators.hmc.DiscontinuousHamiltonianMonteCarloOperatorParser;
import dr.inferencexml.operators.hmc.MixedDiscontinuousHamiltonianMonteCarloOperatorParser;
import dr.xml.*;
import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

import java.io.StringReader;

public class DiscontinuousHMCXMLParserTest extends TestCase {

    public void testLikelihoodBasedPotentialProviderRestoresParameter() {
        Parameter parameter = new Parameter.Default(new double[]{0.2, 0.25});
        ToyMixedLikelihood likelihood = new ToyMixedLikelihood(parameter);
        LikelihoodBasedDiscontinuousPotentialProvider provider =
                new LikelihoodBasedDiscontinuousPotentialProvider(likelihood, parameter);

        double originalX = parameter.getParameterValue(0);
        double originalY = parameter.getParameterValue(1);
        double proposed = provider.getLogDensityAfterSingleCoordinateMove(1, 1.25);

        assertEquals(originalX, parameter.getParameterValue(0), 1E-12);
        assertEquals(originalY, parameter.getParameterValue(1), 1E-12);
        assertFalse(Double.isInfinite(proposed));
    }

    public void testLikelihoodBasedPotentialProviderInvalidatesCacheAcrossRepeatedProbes() {
        Parameter parameter = new Parameter.Default(new double[]{0.0, 0.25});
        ToyMixedLikelihood likelihood = new ToyMixedLikelihood(parameter);
        LikelihoodBasedDiscontinuousPotentialProvider provider =
                new LikelihoodBasedDiscontinuousPotentialProvider(likelihood, parameter);

        double original = provider.getLogDensity();
        assertEquals(Math.log(0.25), original, 1E-12);

        double proposedLow = provider.getLogDensityAfterSingleCoordinateMove(1, 0.5);
        assertEquals(Math.log(0.25), proposedLow, 1E-12);

        double proposedHigh = provider.getLogDensityAfterSingleCoordinateMove(1, 1.5);
        assertEquals(Math.log(0.75), proposedHigh, 1E-12);

        assertEquals(original, provider.getLogDensity(), 1E-12);
    }

    public void testLikelihoodBasedPotentialProviderDefaultUsesParameterEvents() {
        Parameter parameter = new Parameter.Default(new double[]{0.25});
        EventAwareLikelihood likelihood = new EventAwareLikelihood(parameter);
        LikelihoodBasedDiscontinuousPotentialProvider provider =
                new LikelihoodBasedDiscontinuousPotentialProvider(likelihood, parameter);

        assertFalse(provider.isQuiet());
        assertEquals(1.25, provider.getLogDensityAfterSingleCoordinateMove(0, 1.25), 1E-12);
        assertEquals(0.25, parameter.getParameterValue(0), 1E-12);
        assertEquals(0.25, provider.getLogDensity(), 1E-12);
    }

    public void testLikelihoodBasedPotentialProviderQuietSkipsParameterEvents() {
        Parameter parameter = new Parameter.Default(new double[]{0.25});
        EventAwareLikelihood likelihood = new EventAwareLikelihood(parameter);
        LikelihoodBasedDiscontinuousPotentialProvider provider =
                new LikelihoodBasedDiscontinuousPotentialProvider(likelihood, parameter, true);

        assertTrue(provider.isQuiet());
        assertEquals(0.25, provider.getLogDensityAfterSingleCoordinateMove(0, 1.25), 1E-12);
        assertEquals(0.25, parameter.getParameterValue(0), 1E-12);
        assertEquals(0.25, provider.getLogDensity(), 1E-12);
    }

    public void testParseDiscontinuousOperatorXML() throws Exception {
        String xml = "<beast>\n" +
                "  <parameter id=\"p\" value=\"0.25\" lower=\"0.0\" upper=\"2.0\"/>\n" +
                "  <toyTwoStateLikelihood id=\"like\"><parameter idref=\"p\"/></toyTwoStateLikelihood>\n" +
                "  <discontinuousHamiltonianMonteCarloOperator weight=\"1.0\" stepSize=\"1.0\"\n" +
                "                                                   randomStepSizeFraction=\"0.5\" nSteps=\"1\">\n" +
                "    <parameter idref=\"p\"/>\n" +
                "    <likelihoodBasedDiscontinuousPotential quiet=\"true\">\n" +
                "      <toyTwoStateLikelihood idref=\"like\"/>\n" +
                "      <parameter idref=\"p\"/>\n" +
                "    </likelihoodBasedDiscontinuousPotential>\n" +
                "  </discontinuousHamiltonianMonteCarloOperator>\n" +
                "</beast>";

        XMLParser parser = makeParser();
        DiscontinuousHamiltonianMonteCarloOperator parsed =
                (DiscontinuousHamiltonianMonteCarloOperator) parser.parse(new StringReader(xml),
                        DiscontinuousHamiltonianMonteCarloOperator.class);

        assertTrue(parsed instanceof DiscontinuousHamiltonianMonteCarloOperator);
        assertEquals(0.5, parsed.getRandomStepSizeFraction(), 1E-12);
    }

    public void testNegativeRandomStepSizeFractionIsRejected() throws Exception {
        String xml = "<beast>\n" +
                "  <parameter id=\"p\" value=\"0.25\" lower=\"0.0\" upper=\"2.0\"/>\n" +
                "  <toyTwoStateLikelihood id=\"like\"><parameter idref=\"p\"/></toyTwoStateLikelihood>\n" +
                "  <discontinuousHamiltonianMonteCarloOperator weight=\"1.0\" stepSize=\"1.0\"\n" +
                "                                                   randomStepSizeFraction=\"-0.5\" nSteps=\"1\">\n" +
                "    <parameter idref=\"p\"/>\n" +
                "    <likelihoodBasedDiscontinuousPotential quiet=\"true\">\n" +
                "      <toyTwoStateLikelihood idref=\"like\"/>\n" +
                "      <parameter idref=\"p\"/>\n" +
                "    </likelihoodBasedDiscontinuousPotential>\n" +
                "  </discontinuousHamiltonianMonteCarloOperator>\n" +
                "</beast>";

        try {
            makeParser().parse(new StringReader(xml), DiscontinuousHamiltonianMonteCarloOperator.class);
            fail("Expected negative randomStepSizeFraction to be rejected");
        } catch (XMLParseException e) {
            assertTrue(e.getMessage().contains("Random step size fraction"));
        }
    }

    public void testParseMixedOperatorXML() throws Exception {
        String xml = "<beast>\n" +
                "  <parameter id=\"p\" value=\"0.0 0.25\" lower=\"-100.0 0.0\" upper=\"100.0 2.0\"/>\n" +
                "  <parameter id=\"mask\" value=\"0 1\"/>\n" +
                "  <toyMixedLikelihood id=\"like\"><parameter idref=\"p\"/></toyMixedLikelihood>\n" +
                "  <toyMixedGradient id=\"grad\"><parameter idref=\"p\"/></toyMixedGradient>\n" +
                "  <mixedDiscontinuousHamiltonianMonteCarloOperator weight=\"1.0\" stepSize=\"0.25\"\n" +
                "                                                   randomStepSizeFraction=\"0.5\" nSteps=\"4\"\n" +
                "                                                   autoOptimize=\"true\" targetAcceptanceProbability=\"0.8\">\n" +
                "    <parameter idref=\"p\"/>\n" +
                "    <toyMixedGradient idref=\"grad\"/>\n" +
                "    <likelihoodBasedDiscontinuousPotential>\n" +
                "      <toyMixedLikelihood idref=\"like\"/>\n" +
                "      <parameter idref=\"p\"/>\n" +
                "    </likelihoodBasedDiscontinuousPotential>\n" +
                "    <mask><parameter idref=\"mask\"/></mask>\n" +
                "  </mixedDiscontinuousHamiltonianMonteCarloOperator>\n" +
                "</beast>";

        XMLParser parser = makeParser();
        MixedDiscontinuousHamiltonianMonteCarloOperator parsed =
                (MixedDiscontinuousHamiltonianMonteCarloOperator) parser.parse(new StringReader(xml),
                        MixedDiscontinuousHamiltonianMonteCarloOperator.class);

        assertTrue(parsed instanceof MixedDiscontinuousHamiltonianMonteCarloOperator);
        assertEquals(AdaptationMode.ADAPTATION_ON, parsed.getMode());
        assertEquals(0.8, parsed.getTargetAcceptanceProbability(), 1E-12);
        assertEquals(0.25, parsed.getRawParameter(), 1E-12);
        assertEquals(0.5, parsed.getRandomStepSizeFraction(), 1E-12);
    }

    public static Test suite() {
        return new TestSuite(DiscontinuousHMCXMLParserTest.class);
    }

    private XMLParser makeParser() {
        XMLParser parser = new XMLParser(true, true, true, null);
        parser.addXMLObjectParser(new AttributeParser());
        parser.addXMLObjectParser(new ParameterParser());
        parser.addXMLObjectParser(new LikelihoodBasedDiscontinuousPotentialProviderParser());
        parser.addXMLObjectParser(new DiscontinuousHamiltonianMonteCarloOperatorParser());
        parser.addXMLObjectParser(new MixedDiscontinuousHamiltonianMonteCarloOperatorParser());
        parser.addXMLObjectParser(new ToyTwoStateLikelihoodParser());
        parser.addXMLObjectParser(new ToyMixedLikelihoodParser());
        parser.addXMLObjectParser(new ToyMixedGradientParser());
        return parser;
    }

    private static class ToyTwoStateLikelihood extends Likelihood.Abstract {
        private final Parameter parameter;

        ToyTwoStateLikelihood(Parameter parameter) {
            super(null);
            this.parameter = parameter;
        }

        @Override
        protected double calculateLogLikelihood() {
            double y = parameter.getParameterValue(0);
            if (y < 0.0 || y > 2.0) {
                return Double.NEGATIVE_INFINITY;
            }
            return y <= 1.0 ? Math.log(0.25) : Math.log(0.75);
        }
    }

    private static class ToyMixedLikelihood extends Likelihood.Abstract {
        private final Parameter parameter;

        ToyMixedLikelihood(Parameter parameter) {
            super(null);
            this.parameter = parameter;
        }

        @Override
        protected double calculateLogLikelihood() {
            double[] values = parameter.getParameterValues();
            double x = values[0];
            double y = values[1];
            if (y < 0.0 || y > 2.0) {
                return Double.NEGATIVE_INFINITY;
            }
            double discrete = y <= 1.0 ? Math.log(0.25) : Math.log(0.75);
            return -0.5 * x * x + discrete;
        }
    }

    private static class ToyMixedGradient implements GradientWrtParameterProvider {
        private final Parameter parameter;
        private final Likelihood likelihood;

        ToyMixedGradient(Parameter parameter) {
            this.parameter = parameter;
            this.likelihood = new ToyMixedLikelihood(parameter);
        }

        @Override
        public Likelihood getLikelihood() {
            return likelihood;
        }

        @Override
        public Parameter getParameter() {
            return parameter;
        }

        @Override
        public int getDimension() {
            return parameter.getDimension();
        }

        @Override
        public double[] getGradientLogDensity() {
            return new double[]{-parameter.getParameterValue(0), 0.0};
        }
    }

    private static class EventAwareLikelihood extends Likelihood.Abstract implements VariableListener {
        private final Parameter parameter;
        private double cachedValue;

        EventAwareLikelihood(Parameter parameter) {
            super(null);
            this.parameter = parameter;
            this.cachedValue = parameter.getParameterValue(0);
            parameter.addParameterListener(this);
        }

        @Override
        public void variableChangedEvent(Variable variable, int index, Variable.ChangeType type) {
            cachedValue = parameter.getParameterValue(0);
            makeDirty();
        }

        @Override
        protected double calculateLogLikelihood() {
            return cachedValue;
        }
    }

    private static class ToyTwoStateLikelihoodParser extends AbstractXMLObjectParser {
        @Override
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            return new ToyTwoStateLikelihood((Parameter) xo.getChild(Parameter.class));
        }

        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            return new XMLSyntaxRule[]{ new ElementRule(Parameter.class) };
        }

        @Override
        public String getParserDescription() {
            return null;
        }

        @Override
        public Class getReturnType() {
            return Likelihood.class;
        }

        @Override
        public String getParserName() {
            return "toyTwoStateLikelihood";
        }
    }

    private static class ToyMixedLikelihoodParser extends AbstractXMLObjectParser {
        @Override
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            return new ToyMixedLikelihood((Parameter) xo.getChild(Parameter.class));
        }

        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            return new XMLSyntaxRule[]{ new ElementRule(Parameter.class) };
        }

        @Override
        public String getParserDescription() {
            return null;
        }

        @Override
        public Class getReturnType() {
            return Likelihood.class;
        }

        @Override
        public String getParserName() {
            return "toyMixedLikelihood";
        }
    }

    private static class ToyMixedGradientParser extends AbstractXMLObjectParser {
        @Override
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            return new ToyMixedGradient((Parameter) xo.getChild(Parameter.class));
        }

        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            return new XMLSyntaxRule[]{ new ElementRule(Parameter.class) };
        }

        @Override
        public String getParserDescription() {
            return null;
        }

        @Override
        public Class getReturnType() {
            return GradientWrtParameterProvider.class;
        }

        @Override
        public String getParserName() {
            return "toyMixedGradient";
        }
    }
}
