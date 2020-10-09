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

import dr.inference.distribution.shrinkage.JointBayesianBridgeDistributionModel;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.hmc.HessianWrtParameterProvider;
import dr.inference.hmc.ReversibleHMCProvider;
import dr.inference.model.Parameter;
import dr.inference.operators.AdaptableMCMCOperator;
import dr.inference.operators.AdaptationMode;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.hmc.HamiltonianMonteCarloOperator;
import dr.inference.operators.hmc.MassPreconditioner;
import dr.util.Transform;
import dr.xml.*;

import static dr.util.Transform.Util.parseTransform;

/**
 * @author Max Tolkoff
 * @author Marc A. Suchard
 */

public class HamiltonianMonteCarloOperatorParser extends AbstractXMLObjectParser {

    private final static String HMC_OPERATOR = "hamiltonianMonteCarloOperator";
    private final static String N_STEPS = "nSteps";
    private final static String STEP_SIZE = "stepSize";
    private final static String RANDOM_STEP_FRACTION = "randomStepCountFraction";
    private final static String PRECONDITIONING = "preconditioning";
    private final static String PRECONDITIONING_UPDATE_FREQUENCY = "preconditioningUpdateFrequency";
    private final static String PRECONDITIONING_DELAY = "preconditioningDelay";
    private final static String PRECONDITIONING_MEMORY = "preconditioningMemory";
    private final static String PRECONDITIONER = "preconditioner";
    private final static String SHRINKAGE_PRECONDITIONER = "shrinkagePreconditioner";
    private final static String GRADIENT_CHECK_COUNT = "gradientCheckCount";
    public final static String GRADIENT_CHECK_TOLERANCE = "gradientCheckTolerance";
    private final static String MAX_ITERATIONS = "checkStepSizeMaxIterations";
    private final static String REDUCTION_FACTOR = "checkStepSizeReductionFactor";
    private final static String TARGET_ACCEPTANCE_PROBABILITY = "targetAcceptanceProbability";
    private final static String INSTABILITY_HANDLER = "instabilityHandler";
    private final static String MASK = "mask";

    private MassPreconditioner preconditioner;

    @Override
    public String getParserName() {
        return HMC_OPERATOR;
    }

    private MassPreconditioner.Type parsePreconditioning(XMLObject xo) throws XMLParseException {

        return MassPreconditioner.Type.parseFromString(
                xo.getAttribute(PRECONDITIONING, MassPreconditioner.Type.NONE.getName())
        );
    }

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);
        int nSteps = xo.getAttribute(N_STEPS, 10);
        double stepSize = xo.getDoubleAttribute(STEP_SIZE);

        MassPreconditioner.Type preconditioningType = parsePreconditioning(xo);

        double randomStepFraction = Math.abs(xo.getAttribute(RANDOM_STEP_FRACTION, 0.0));
        if (randomStepFraction > 1) {
            throw new XMLParseException("Random step count fraction must be < 1.0");
        }

        int preconditioningUpdateFrequency = xo.getAttribute(PRECONDITIONING_UPDATE_FREQUENCY, 0);

        int preconditioningDelay = xo.getAttribute(PRECONDITIONING_DELAY, 0);

        int preconditioningMemory = xo.getAttribute(PRECONDITIONING_MEMORY, 0);

        AdaptationMode adaptationMode = AdaptationMode.parseMode(xo);

        GradientWrtParameterProvider derivative =
                (GradientWrtParameterProvider) xo.getChild(GradientWrtParameterProvider.class);

        if (preconditioningType != MassPreconditioner.Type.NONE &&
                !(derivative instanceof HessianWrtParameterProvider)) {
            throw new XMLParseException("Unable precondition without a Hessian provider");
        }

        Parameter parameter = (Parameter) xo.getChild(Parameter.class);

        if (parameter == null) {
            parameter = derivative.getParameter();
        }

        Transform transform = parseTransform(xo);

        ReversibleHMCProvider reversibleHMCprovider = (ReversibleHMCProvider) xo.getChild(ReversibleHMCProvider.class);

        boolean dimensionMismatch = derivative.getDimension() != parameter.getDimension();
        if (transform instanceof Transform.MultivariableTransform) {
            dimensionMismatch = ((Transform.MultivariableTransform) transform).getDimension() != parameter.getDimension();
        }

        if (dimensionMismatch) {
            throw new XMLParseException("Gradient (" + derivative.getDimension() +
                    ") must be the same dimensions as the parameter (" + parameter.getDimension() + ")");
        }

        Parameter mask = null;
        if (xo.hasChildNamed(MASK)) {
            mask = (Parameter) xo.getElementFirstChild(MASK);

            if (mask.getDimension() != derivative.getDimension()) {
                throw new XMLParseException("Mask (" + mask.getDimension()
                        + ") must be the same dimension as the gradient (" + derivative.getDimension() + ")");
            }
        }

        int gradientCheckCount = xo.getAttribute(GRADIENT_CHECK_COUNT, 0);
        double gradientCheckTolerance = xo.getAttribute(GRADIENT_CHECK_TOLERANCE, 1E-3);
        int maxIterations = xo.getAttribute(MAX_ITERATIONS, 10);
        double reductionFactor = xo.getAttribute(REDUCTION_FACTOR, 0.1);
        double targetAcceptanceProbability = xo.getAttribute(TARGET_ACCEPTANCE_PROBABILITY,
                0.8); // Stan default
        String instabilityHandlerCase = xo.getAttribute(INSTABILITY_HANDLER, "reject");
        HamiltonianMonteCarloOperator.InstabilityHandler instabilityHandler = HamiltonianMonteCarloOperator.InstabilityHandler.factory(instabilityHandlerCase);

        HamiltonianMonteCarloOperator.Options runtimeOptions = new HamiltonianMonteCarloOperator.Options(
                stepSize, nSteps, randomStepFraction,
                preconditioningUpdateFrequency, preconditioningDelay, preconditioningMemory,
                gradientCheckCount, gradientCheckTolerance,
                maxIterations, reductionFactor,
                targetAcceptanceProbability,
                instabilityHandler
        );

        if (xo.hasChildNamed(PRECONDITIONER)) {
            if (preconditioningType != MassPreconditioner.Type.NONE) {
                throw new XMLParseException("Cannot precondition and use an alternative preconditioner");
            }
            XMLObject cxo = xo.getChild(PRECONDITIONER);

            if (cxo.hasChildNamed(SHRINKAGE_PRECONDITIONER)) {
                XMLObject ccxo = cxo.getChild(SHRINKAGE_PRECONDITIONER);
                JointBayesianBridgeDistributionModel bridge = (JointBayesianBridgeDistributionModel) ccxo.getChild(JointBayesianBridgeDistributionModel.class);
                preconditioningUpdateFrequency = 1;
                preconditioner = new MassPreconditioner.ShrinkagePreconditioner(bridge, transform);
            } else {
                throw new XMLParseException("Unknown preconditioner specified");
            }
        } else {
            preconditioner = preconditioningType.factory(derivative, transform, runtimeOptions);
        }


        if (xo.hasChildNamed(PRECONDITIONER)) {
            return new HamiltonianMonteCarloOperator(adaptationMode, weight, derivative,
                    parameter, transform, mask,
                    runtimeOptions, preconditioner);
        } else {
            return factory(adaptationMode, weight, derivative, parameter, transform, mask, runtimeOptions, preconditioner, reversibleHMCprovider);
        }

    }

    protected HamiltonianMonteCarloOperator factory(AdaptationMode adaptationMode, double weight, GradientWrtParameterProvider derivative,
                                                    Parameter parameter, Transform transform, Parameter mask,
                                                    HamiltonianMonteCarloOperator.Options runtimeOptions, MassPreconditioner preconditioner,
                                                    ReversibleHMCProvider reversibleHMCprovider) {

        return new HamiltonianMonteCarloOperator(adaptationMode, weight, derivative,
                parameter, transform, mask,
                runtimeOptions, preconditioner);

    }


    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    protected final XMLSyntaxRule[] rules = {
            AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
            AttributeRule.newIntegerRule(N_STEPS, true),
            AttributeRule.newDoubleRule(STEP_SIZE),
            AttributeRule.newBooleanRule(AdaptableMCMCOperator.AUTO_OPTIMIZE, true),
            AttributeRule.newStringRule(PRECONDITIONING, true),
            AttributeRule.newDoubleRule(RANDOM_STEP_FRACTION, true),
            AttributeRule.newDoubleRule(TARGET_ACCEPTANCE_PROBABILITY, true),
            new ElementRule(Parameter.class, true),
            new ElementRule(Transform.MultivariableTransformWithParameter.class, true),
            new ElementRule(GradientWrtParameterProvider.class),
            new ElementRule(MASK, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class),

            }, true),
            new ElementRule(PRECONDITIONER, new XMLSyntaxRule[]{
                    new ElementRule(MassPreconditioner.class)
            }, true),

    };

    @Override
    public String getParserDescription() {
        return "Returns a Hamiltonian Monte Carlo transition kernel";
    }

    @Override
    public Class getReturnType() {
        return HamiltonianMonteCarloOperator.class;
    }
}
