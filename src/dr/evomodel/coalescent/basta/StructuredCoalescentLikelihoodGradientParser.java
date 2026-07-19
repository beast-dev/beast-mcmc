/*
 * StructuredCoalescentLikelihoodGradientParser.java
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

package dr.evomodel.coalescent.basta;

import dr.evomodel.coalescent.AbstractStructuredCoalescentLikelihood;
import dr.evomodel.coalescent.mascot.MascotGradient;
import dr.evomodel.coalescent.mascot.MascotLikelihood;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.hmc.NumericGradientStepSizeProvider;
import dr.xml.*;

/**
 * Parses {@code <structuredCoalescentLikelihoodGradient>}, dispatching on
 * the Java type of the referenced likelihood (a {@link BastaLikelihood} or a
 * {@link MascotLikelihood}) rather than a separate {@code type} attribute --
 * unlike {@link StructuredCoalescentLikelihoodParser}, the likelihood object
 * is already resolved by the time this parser runs, so there's nothing to
 * dispatch on except that.
 * <p/>
 * {@code wrtParameter} keeps BASTA's existing vocabulary
 * ({@code "migrationRate"}/{@code "populationSize"}) for both engines;
 * MASCOT additionally accepts {@code "clockRate"}, a capability BASTA
 * doesn't have (its likelihood currently requires a {@code
 * StrictClockBranchRates} clock with no differentiable gradient). BASTA's
 * required {@code SubstitutionModel} child (used for its migration-rate
 * chain rule) is not required for MASCOT, which doesn't need one.
 */
public class StructuredCoalescentLikelihoodGradientParser extends AbstractXMLObjectParser {

    private static final String NAME = "structuredCoalescentLikelihoodGradient";
    private static final String WRT_PARAMETER = "wrtParameter";
    private static final String MIGRATION_RATE = "migrationRate";
    private static final String POPULATION_SIZE = "populationSize";
    private static final String CLOCK_RATE = "clockRate";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        AbstractStructuredCoalescentLikelihood likelihood =
                (AbstractStructuredCoalescentLikelihood) xo.getChild(AbstractStructuredCoalescentLikelihood.class);
        String wrtParameter = xo.getStringAttribute(WRT_PARAMETER);

        if (likelihood instanceof MascotLikelihood) {
            return parseMascotGradient((MascotLikelihood) likelihood, wrtParameter, xo);
        }
        return parseBastaGradient((BastaLikelihood) likelihood, wrtParameter, xo);
    }

    private Object parseBastaGradient(BastaLikelihood likelihood, String wrtParameter, XMLObject xo) throws XMLParseException {
        SubstitutionModel substitutionModel = (SubstitutionModel) xo.getChild(SubstitutionModel.class);
        if (substitutionModel == null) {
            throw new XMLParseException(NAME + " requires a SubstitutionModel child for a BASTA likelihood");
        }
        StructuredCoalescentLikelihoodGradient.WrtParameter type =
                StructuredCoalescentLikelihoodGradient.WrtParameter.factory(wrtParameter);
        if (type == null) {
            throw new XMLParseException(WRT_PARAMETER + " must be \"" + MIGRATION_RATE + "\" or \"" +
                    POPULATION_SIZE + "\" for a BASTA likelihood, found \"" + wrtParameter + "\"");
        }
        return new StructuredCoalescentLikelihoodGradient(likelihood, substitutionModel, type);
    }

    private Object parseMascotGradient(MascotLikelihood likelihood, String wrtParameter, XMLObject xo) throws XMLParseException {
        MascotGradient.Part part;
        if (wrtParameter.equalsIgnoreCase(MIGRATION_RATE)) {
            part = MascotGradient.Part.MIGRATION;
        } else if (wrtParameter.equalsIgnoreCase(POPULATION_SIZE)) {
            part = MascotGradient.Part.POPULATION_SIZE;
        } else if (wrtParameter.equalsIgnoreCase(CLOCK_RATE)) {
            part = MascotGradient.Part.CLOCK_RATE;
        } else {
            throw new XMLParseException(WRT_PARAMETER + " must be \"" + MIGRATION_RATE + "\", \"" +
                    POPULATION_SIZE + "\", or \"" + CLOCK_RATE + "\" for a MASCOT likelihood, found \"" +
                    wrtParameter + "\"");
        }
        MascotGradient gradient = new MascotGradient(likelihood, part);
        gradient.setNumericGradientStepSize(NumericGradientStepSizeProvider.parseStepSizeRatio(xo));
        return gradient;
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newStringRule(WRT_PARAMETER),
            new ElementRule(AbstractStructuredCoalescentLikelihood.class),
            // Required for a BASTA likelihood, unused for a MASCOT one --
            // declared optional here and enforced imperatively in
            // parseBastaGradient, the same pattern used for BASTA's
            // popSizes/logPopSizes mutual exclusivity in
            // StructuredCoalescentLikelihoodParser.
            new ElementRule(SubstitutionModel.class, true),
    };

    @Override
    public String getParserDescription() {
        return "Gradient of a structured-coalescent likelihood (BASTA or MASCOT) with respect to its " +
                "migration-rate or population-size parameter (wrtParameter=\"migrationRate\" or " +
                "\"populationSize\"), or, for a MASCOT likelihood only, its clock rate " +
                "(wrtParameter=\"clockRate\").";
    }

    @Override
    public Class getReturnType() {
        return GradientWrtParameterProvider.class;
    }

    @Override
    public String getParserName() {
        return NAME;
    }
}
