/*
 * HamiltonianMonteCarloOperatorParser.java
 *
 * Copyright (c) 2002-2017 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.inferencexml.operators.hmc;

import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.hmc.PrecisionMatrixVectorProductProvider;
import dr.inference.model.Parameter;
import dr.inference.operators.CoercableMCMCOperator;
import dr.inference.operators.CoercionMode;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.hmc.BouncyParticleOperator;
import dr.xml.*;

/**
 * @author Zhenyu Zhang
 * @author Marc A. Suchard
 */

public class BouncyParticleOperatorParser extends AbstractXMLObjectParser {

    private final static String BPO_OPERATOR = "bouncyParticleOperator";
    private final static String RANDOM_TIME_WIDTH = "randomTimeWidth";
    private final static String UPDATE_FREQUENCY = "preconditioningUpdateFrequency";
    private final static String TIME_GUESS = "roughTravelTimeGuess";
    private final static String MASKING = "mask";

    @Override
    public String getParserName() {
        return BPO_OPERATOR;
    }

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);

        @SuppressWarnings("unused") CoercionMode coercionMode = CoercionMode.parseMode(xo);

        GradientWrtParameterProvider derivative =
                (GradientWrtParameterProvider) xo.getChild(GradientWrtParameterProvider.class);

        PrecisionMatrixVectorProductProvider productProvider = (PrecisionMatrixVectorProductProvider)
                xo.getChild(PrecisionMatrixVectorProductProvider.class);

        Parameter mask = null;
        if (xo.hasChildNamed(MASKING)) {
            mask = (Parameter) xo.getElementFirstChild(MASKING);
        }

        BouncyParticleOperator.Options runtimeOptions = parseRuntimeOptions(xo);

        return new BouncyParticleOperator(derivative, productProvider, weight, runtimeOptions, mask);
    }

    private BouncyParticleOperator.Options parseRuntimeOptions(XMLObject xo) throws XMLParseException {

        double randomTimeWidth = xo.getAttribute(RANDOM_TIME_WIDTH, 0.5);
        int updateFrequency = xo.getAttribute(UPDATE_FREQUENCY, 0);
        boolean runTimeGuess = xo.getAttribute(TIME_GUESS, false);

        return new BouncyParticleOperator.Options(randomTimeWidth, updateFrequency, runTimeGuess);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
            AttributeRule.newBooleanRule(CoercableMCMCOperator.AUTO_OPTIMIZE, true),
            AttributeRule.newDoubleRule(RANDOM_TIME_WIDTH, true),
            AttributeRule.newIntegerRule(UPDATE_FREQUENCY, true),
            new ElementRule(GradientWrtParameterProvider.class),
            new ElementRule(PrecisionMatrixVectorProductProvider.class),
            new ElementRule(MASKING, new XMLSyntaxRule[] {
                    new ElementRule(Parameter.class),
            }, true),
    };

    @Override
    public String getParserDescription() {
        return "Returns a bouncy particle transition kernel for truncated normals";
    }

    @Override
    public Class getReturnType() {
        return BouncyParticleOperator.class;
    }
}