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
import dr.inference.hmc.PrecisionColumnProvider;
import dr.inference.hmc.PrecisionMatrixVectorProductProvider;
import dr.inference.model.Parameter;
import dr.inference.operators.AdaptableMCMCOperator;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.hmc.AbstractParticleOperator;
import dr.inference.operators.hmc.BouncyParticleOperator;
import dr.inference.operators.hmc.MassPreconditionScheduler;
import dr.xml.*;

/**
 * @author Zhenyu Zhang
 * @author Marc A. Suchard
 */

public class BouncyParticleOperatorParser extends AbstractXMLObjectParser {

    private final static String BPO_OPERATOR = "bouncyParticleOperator";
    private final static String RANDOM_TIME_WIDTH = "randomTimeWidth";
    private final static String UPDATE_FREQUENCY = "preconditioningUpdateFrequency";

    private final static String UPDATE_SCM_DELAY = "updateSampleCovDelay";
    private final static String UPDATE_SCM_FREQUENCY = "updateSampleCovFrequency";

    private final static String MASKING = "mask";
    private final static String REFRESH_VELOCITY = "refreshVelocity";

    private final static String TEST_NATIVE_BOUNCE = "testNativeFindBounce";
    private final static String USE_NATIVE_BOUNCE = "useNativeFindBounce";
    private final static String USE_NATIVE_UPDATE_DYNAMICS = "useNativeUpdateDynamics";

    @Override
    public String getParserName() {
        return BPO_OPERATOR;
    }

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        boolean refreshVelocity = xo.getAttribute(REFRESH_VELOCITY, true);

        return new BouncyParticleOperator(
                (GradientWrtParameterProvider) xo.getChild(GradientWrtParameterProvider.class),
                (PrecisionMatrixVectorProductProvider) xo.getChild(PrecisionMatrixVectorProductProvider.class),
                (PrecisionColumnProvider) xo.getChild(PrecisionColumnProvider.class),
                xo.getDoubleAttribute(MCMCOperator.WEIGHT),
                parseRuntimeOptions(xo),
                parseNativeCodeOptions(xo),
                refreshVelocity,
                parseMask(xo),
                null, MassPreconditionScheduler.Type.NONE);
    }

    static Parameter parseMask(XMLObject xo) throws XMLParseException {
        Parameter mask = null;
        if (xo.hasChildNamed(MASKING)) {
            mask = (Parameter) xo.getElementFirstChild(MASKING);
        }

        return mask;
    }


    static AbstractParticleOperator.Options parseRuntimeOptions(XMLObject xo) throws XMLParseException {

        double randomTimeWidth = xo.getAttribute(RANDOM_TIME_WIDTH, 0.5);
        int updateFrequency = xo.getAttribute(UPDATE_FREQUENCY, 0);
        int preconditioningMaxUpdate = xo.getAttribute(PreconditionHandlerParser.PRECONDITIONING_MAX_UPDATE, 0);
        int preconditioningDelay = xo.getAttribute(PreconditionHandlerParser.PRECONDITIONING_DELAY, 0);

        int updateSampleCovFrequency = xo.getAttribute(UPDATE_SCM_FREQUENCY, 0);
        int updateSampleCovDelay = xo.getAttribute(UPDATE_SCM_DELAY, 0);

        return new AbstractParticleOperator.Options(randomTimeWidth, updateFrequency, preconditioningMaxUpdate,
                preconditioningDelay, updateSampleCovFrequency, updateSampleCovDelay);
    }

    static AbstractParticleOperator.NativeCodeOptions parseNativeCodeOptions(XMLObject xo) throws XMLParseException {

        boolean testNativeFindNextBounce = xo.getAttribute(TEST_NATIVE_BOUNCE, false);
        boolean useNativeFindNextBounce = xo.getAttribute(USE_NATIVE_BOUNCE, false);
        boolean useNativeUpdateDynamics = xo.getAttribute(USE_NATIVE_UPDATE_DYNAMICS, false);

        return new AbstractParticleOperator.NativeCodeOptions(testNativeFindNextBounce, useNativeFindNextBounce, useNativeUpdateDynamics);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    final static XMLSyntaxRule[] rules = {
            AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
            AttributeRule.newBooleanRule(AdaptableMCMCOperator.AUTO_OPTIMIZE, true),
            AttributeRule.newDoubleRule(RANDOM_TIME_WIDTH, true),
            AttributeRule.newIntegerRule(UPDATE_FREQUENCY, true),
            AttributeRule.newIntegerRule(UPDATE_SCM_FREQUENCY, true),
            AttributeRule.newIntegerRule(UPDATE_SCM_DELAY, true),
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