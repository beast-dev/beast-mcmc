/*
 * SplitHamiltonianMonteCarloOperatorParser.java
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

package dr.inferencexml.operators.hmc;

import dr.inference.hmc.ReversibleHMCProvider;
import dr.inference.model.Parameter;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.hmc.SplitHMCtravelTimeMultiplier;
import dr.inference.operators.hmc.SplitHamiltonianMonteCarloOperator;
import dr.xml.*;

/**
 * @author Zhenyu Zhang
 */

public class SplitHamiltonianMonteCarloOperatorParser extends AbstractXMLObjectParser { //todo: merge with HMC parser?

    private final static String SPLIT_HMC = "splitHamiltonianMonteCarloOperator";
    private final static String N_STEPS = "nSteps";
    private final static String N_INNER_STEPS = "nInnerSteps";
    private final static String STEP_SIZE = "stepSize";
    private final static String RELATIVE_SCALE = "relativeScale";
    private final static String UPDATE_RS_FREQUENCY = "updateRelativeScaleFrequency";
    private final static String UPDATE_RS_DELAY = "updateRelativeScaleDelay";
    private final static String GET_RS_DELAY = "getRelativeScaleDelay";
    private final static String GET_RS_FREQUENCY = "getRelativeScaleFrequency";
    private final static String GRADIENT_CHECK_COUNT = "gradientCheckCount";
    private final static String GRADIENT_CHECK_TOL = "gradientCheckTolerance";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);
        double stepSize = xo.getDoubleAttribute(STEP_SIZE);
        double relativeScale = xo.getDoubleAttribute(RELATIVE_SCALE);

        ReversibleHMCProvider reversibleHMCproviderInner = (ReversibleHMCProvider) xo.getChild(1); //todo: avoid hard-coded order of reversible provider?
        ReversibleHMCProvider reversibleHMCproviderOuter = (ReversibleHMCProvider) xo.getChild(2);

        int nStep = xo.getAttribute(N_STEPS, 5);
        int nInnerStep = xo.getAttribute(N_INNER_STEPS, 5);
        int gradientCheckCount = xo.getAttribute(GRADIENT_CHECK_COUNT, 0);
        double gradientCheckTol = xo.getAttribute(GRADIENT_CHECK_TOL, 0.01);

        SplitHMCtravelTimeMultiplier.RSoptions rsOptions = parseRSoptions(xo);
        SplitHMCtravelTimeMultiplier splitHMCmultiplier = SplitHMCtravelTimeMultiplier.create(reversibleHMCproviderInner, reversibleHMCproviderOuter, rsOptions);

        Parameter parameter = (Parameter) xo.getChild(Parameter.class);

        return new SplitHamiltonianMonteCarloOperator(weight, reversibleHMCproviderInner, reversibleHMCproviderOuter, parameter,
                stepSize, relativeScale, nStep
                , nInnerStep, gradientCheckCount, gradientCheckTol, splitHMCmultiplier);
    }

    static SplitHMCtravelTimeMultiplier.RSoptions parseRSoptions(XMLObject xo) throws XMLParseException {

        int updateRSdelay = xo.getAttribute(UPDATE_RS_DELAY, 0);
        int updateRSfrequency = xo.getAttribute(UPDATE_RS_FREQUENCY, 0);
        int getRSdelay = xo.getAttribute(GET_RS_DELAY, 0);
        int getRSfrequency = xo.getAttribute(GET_RS_FREQUENCY, 0);

        return new SplitHMCtravelTimeMultiplier.RSoptions(updateRSdelay, updateRSfrequency, getRSdelay, getRSfrequency);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    protected final XMLSyntaxRule[] rules = {
            AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
            AttributeRule.newIntegerRule(N_STEPS, true),
            AttributeRule.newIntegerRule(N_INNER_STEPS, true),
            AttributeRule.newDoubleRule(STEP_SIZE),
            AttributeRule.newDoubleRule(RELATIVE_SCALE)
    };

    @Override
    public String getParserDescription() {
        return null;
    }

    @Override
    public Class getReturnType() {
        return null;
    }

    @Override
    public String getParserName() {
        return SPLIT_HMC;
    }
}
