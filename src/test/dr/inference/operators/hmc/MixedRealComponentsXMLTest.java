/*
 * MixedRealComponentsXMLTest.java
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

import dr.inference.distribution.EmbeddedOrdinalLikelihood;
import dr.inference.hmc.CompoundDiscontinuousPotentialProvider;
import dr.inference.hmc.CompoundGradient;
import dr.inference.model.CompoundLikelihood;
import dr.inference.model.ParameterParser;
import dr.inference.operators.hmc.MixedDiscontinuousHamiltonianMonteCarloOperator;
import dr.inferencexml.distribution.DistributionLikelihoodParser;
import dr.inferencexml.distribution.EmbeddedOrdinalLikelihoodParser;
import dr.inferencexml.distribution.NormalDistributionModelParser;
import dr.inferencexml.hmc.CompoundDiscontinuousPotentialProviderParser;
import dr.inferencexml.hmc.CompoundGradientParser;
import dr.inferencexml.hmc.GradientWrapperParser;
import dr.inferencexml.hmc.LikelihoodBasedDiscontinuousPotentialProviderParser;
import dr.inferencexml.hmc.ZeroGradientParser;
import dr.inferencexml.model.CompoundLikelihoodParser;
import dr.inferencexml.operators.hmc.MixedDiscontinuousHamiltonianMonteCarloOperatorParser;
import dr.xml.AttributeParser;
import dr.xml.XMLParser;
import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

import java.io.StringReader;

public class MixedRealComponentsXMLTest extends TestCase {

    public void testParseMixedOperatorFromRealComponentsXML() throws Exception {
        XMLParser parser = makeParser();
        Object parsed = parser.parse(new StringReader(makeXML()), MixedDiscontinuousHamiltonianMonteCarloOperator.class);
        assertTrue(parsed instanceof MixedDiscontinuousHamiltonianMonteCarloOperator);
    }

    public void testRealComponentObjectsComposeCorrectly() throws Exception {
        XMLParser parser = makeParser();
        parser.parse(new StringReader(makeXML()), false);

        CompoundGradient gradient = (CompoundGradient) parser.getObjectStore().get("jointGrad").getNativeObject();
        CompoundDiscontinuousPotentialProvider provider =
                (CompoundDiscontinuousPotentialProvider) parser.getObjectStore().get("jointDisc").getNativeObject();
        CompoundLikelihood likelihood = (CompoundLikelihood) parser.getObjectStore().get("jointLike").getNativeObject();
        EmbeddedOrdinalLikelihood embedded =
                (EmbeddedOrdinalLikelihood) parser.getObjectStore().get("embeddedLike").getNativeObject();

        assertEquals(2, gradient.getDimension());
        assertEquals(2, provider.getDimension());
        assertEquals(0, embedded.getOrdinalState());
        assertEquals(likelihood.getLogLikelihood(), provider.getLogDensity(), 1E-12);
    }

    public static Test suite() {
        return new TestSuite(MixedRealComponentsXMLTest.class);
    }

    private XMLParser makeParser() {
        XMLParser parser = new XMLParser(true, true, true, null);
        parser.addXMLObjectParser(new AttributeParser());
        parser.addXMLObjectParser(new ParameterParser());
        parser.addXMLObjectParser(new NormalDistributionModelParser());
        parser.addXMLObjectParser(new DistributionLikelihoodParser());
        parser.addXMLObjectParser(new GradientWrapperParser());
        parser.addXMLObjectParser(new ZeroGradientParser());
        parser.addXMLObjectParser(new CompoundGradientParser());
        parser.addXMLObjectParser(new EmbeddedOrdinalLikelihoodParser());
        parser.addXMLObjectParser(new CompoundLikelihoodParser());
        parser.addXMLObjectParser(new LikelihoodBasedDiscontinuousPotentialProviderParser());
        parser.addXMLObjectParser(new CompoundDiscontinuousPotentialProviderParser());
        parser.addXMLObjectParser(new MixedDiscontinuousHamiltonianMonteCarloOperatorParser());
        return parser;
    }

    private String makeXML() {
        return "<beast>\n" +
                "  <parameter id=\"x\" value=\"0.0\"/>\n" +
                "  <parameter id=\"latent\" value=\"0.25\" lower=\"0.0\" upper=\"3.0\"/>\n" +
                "  <parameter id=\"zeroMean\" value=\"0.0\"/>\n" +
                "  <parameter id=\"unitStdev\" value=\"1.0\"/>\n" +
                "  <parameter id=\"mask\" value=\"0 1\"/>\n" +
                "  <parameter id=\"logWeights\" value=\"" + Math.log(0.2) + " " + Math.log(0.3) + " " + Math.log(0.5) + "\"/>\n" +
                "  <parameter id=\"cuts\" value=\"0.0 1.0 2.0 3.0\"/>\n" +
                "\n" +
                "  <normalDistributionModel id=\"normalModel\">\n" +
                "    <mean><parameter idref=\"zeroMean\"/></mean>\n" +
                "    <stdev><parameter idref=\"unitStdev\"/></stdev>\n" +
                "  </normalDistributionModel>\n" +
                "\n" +
                "  <distributionLikelihood id=\"xLike\">\n" +
                "    <distribution><normalDistributionModel idref=\"normalModel\"/></distribution>\n" +
                "    <data><parameter idref=\"x\"/></data>\n" +
                "  </distributionLikelihood>\n" +
                "\n" +
                "  <gradient id=\"xGrad\">\n" +
                "    <distributionLikelihood idref=\"xLike\"/>\n" +
                "    <parameter idref=\"x\"/>\n" +
                "  </gradient>\n" +
                "\n" +
                "  <embeddedOrdinalLikelihood id=\"embeddedLike\">\n" +
                "    <latent><parameter idref=\"latent\"/></latent>\n" +
                "    <logWeights><parameter idref=\"logWeights\"/></logWeights>\n" +
                "    <cuts><parameter idref=\"cuts\"/></cuts>\n" +
                "  </embeddedOrdinalLikelihood>\n" +
                "\n" +
                "  <zeroGradient id=\"latentZeroGrad\">\n" +
                "    <parameter idref=\"latent\"/>\n" +
                "    <embeddedOrdinalLikelihood idref=\"embeddedLike\"/>\n" +
                "  </zeroGradient>\n" +
                "\n" +
                "  <compoundGradient id=\"jointGrad\">\n" +
                "    <gradient idref=\"xGrad\"/>\n" +
                "    <zeroGradient idref=\"latentZeroGrad\"/>\n" +
                "  </compoundGradient>\n" +
                "\n" +
                "  <compoundLikelihood id=\"jointLike\">\n" +
                "    <distributionLikelihood idref=\"xLike\"/>\n" +
                "    <embeddedOrdinalLikelihood idref=\"embeddedLike\"/>\n" +
                "  </compoundLikelihood>\n" +
                "\n" +
                "  <likelihoodBasedDiscontinuousPotential id=\"xDiscLike\">\n" +
                "    <distributionLikelihood idref=\"xLike\"/>\n" +
                "    <parameter idref=\"x\"/>\n" +
                "  </likelihoodBasedDiscontinuousPotential>\n" +
                "\n" +
                "  <compoundDiscontinuousPotential id=\"jointDisc\">\n" +
                "    <likelihoodBasedDiscontinuousPotential idref=\"xDiscLike\"/>\n" +
                "    <embeddedOrdinalLikelihood idref=\"embeddedLike\"/>\n" +
                "  </compoundDiscontinuousPotential>\n" +
                "\n" +
                "  <mixedDiscontinuousHamiltonianMonteCarloOperator weight=\"1.0\" stepSize=\"0.25\" nSteps=\"4\">\n" +
                "    <compoundGradient idref=\"jointGrad\"/>\n" +
                "    <compoundDiscontinuousPotential idref=\"jointDisc\"/>\n" +
                "    <mask><parameter idref=\"mask\"/></mask>\n" +
                "  </mixedDiscontinuousHamiltonianMonteCarloOperator>\n" +
                "</beast>";
    }
}
