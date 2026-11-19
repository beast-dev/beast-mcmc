/*
 * CollapsedBridgeHamiltonianOperatorParser.java
 *
 * Copyright © 2002-2026 the BEAST Development Team
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
 * BEAST is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package dr.inferencexml.operators.hmc;

import dr.inference.distribution.shrinkage.BayesianBridgeLikelihood;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.hmc.JointGradient;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.hmc.CollapsedBridgeHamiltonianOperator;
import dr.xml.*;

import java.util.Collections;
import java.util.Set;

/**
 * @author Marc A. Suchard
 */
public class CollapsedBridgeHamiltonianOperatorParser extends AbstractXMLObjectParser {

    private static final String N_STEPS = "nSteps";
    private static final String STEP_SIZE = "stepSize";
    private static final String KAPPA = "kappa";
    private static final String SCORE_TERMS = "scoreTerms";
    private static final String MASS_VARIANCE = "massVariance";

    public String getParserName() {
        return CollapsedBridgeHamiltonianOperator.OPERATOR_NAME;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);
        int nSteps = xo.getIntegerAttribute(N_STEPS);
        double stepSize = xo.getDoubleAttribute(STEP_SIZE);
        double kappa = xo.getAttribute(KAPPA, 10.0);
        if (nSteps < 1 || !(stepSize > 0.0) || !(kappa >= 0.0)) {
            throw new XMLParseException("Need nSteps >= 1, stepSize > 0 and kappa >= 0");
        }

        Parameter rates = (Parameter) xo.getChild(Parameter.class);

        BayesianBridgeLikelihood bridge = null;
        GradientWrtParameterProvider likelihoodGradient = null;
        for (int i = 0; i < xo.getChildCount(); i++) {
            Object child = xo.getChild(i);
            if (child instanceof BayesianBridgeLikelihood) {
                bridge = (BayesianBridgeLikelihood) child;
            } else if (child instanceof GradientWrtParameterProvider) {
                likelihoodGradient = (GradientWrtParameterProvider) child;
            }
        }

        if (bridge.getParameter() != rates) {
            throw new XMLParseException("The bridge must be on parameter " + rates.getId());
        }
        if (bridge.getLocalScale() == null) {
            throw new XMLParseException("The bridge needs local scales (a joint Bayesian bridge)");
        }
        if (likelihoodGradient != null) {
            if (likelihoodGradient instanceof JointGradient) {
                throw new XMLParseException("Give the likelihood-only gradient; the operator adds the bridge terms");
            }
            if (likelihoodGradient.getParameter() != rates
                    || likelihoodGradient.getDimension() != rates.getDimension()) {
                throw new XMLParseException("The gradient must be with respect to " + rates.getId());
            }
        }

        Likelihood scoreTerms = (Likelihood) xo.getElementFirstChild(SCORE_TERMS);
        Set<Likelihood> terms = scoreTerms.getLikelihoodSet();
        if (!terms.contains(bridge)) {
            throw new XMLParseException("scoreTerms must contain the bridge " + bridge.getId()
                    + " (give the prior: the joint minus the data likelihood)");
        }
        if (likelihoodGradient != null && likelihoodGradient.getLikelihood() != null
                && !Collections.disjoint(terms, likelihoodGradient.getLikelihood().getLikelihoodSet())) {
            throw new XMLParseException("scoreTerms must not contain the data likelihood of the gradient "
                    + "(give the prior: the joint minus the data likelihood)");
        }

        Parameter variance = (Parameter) xo.getElementFirstChild(MASS_VARIANCE);
        if (variance.getDimension() != rates.getDimension()) {
            throw new XMLParseException("massVariance has dimension " + variance.getDimension()
                    + ", not " + rates.getDimension());
        }
        double[] massVariance = variance.getParameterValues();
        for (double v : massVariance) {
            if (!(v > 0.0 && v < Double.POSITIVE_INFINITY)) {
                throw new XMLParseException("massVariance values must be positive and finite");
            }
        }

        return new CollapsedBridgeHamiltonianOperator(rates, bridge, likelihoodGradient, scoreTerms,
                massVariance, nSteps, stepSize, kappa, weight);
    }

    public String getParserDescription() {
        return "HMC on Bayesian-bridge coefficients with the local scales integrated out, in warped log " +
                "coordinates, followed by an exact tilted-stable redraw of every local scale. scoreTerms must be " +
                "the joint minus the data likelihood (the prior). Valid only at chain temperature 1: not for " +
                "heated MC3 chains, where the subtracted change in scoreTerms is not tempered";
    }

    public Class getReturnType() {
        return CollapsedBridgeHamiltonianOperator.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
            AttributeRule.newIntegerRule(N_STEPS),
            AttributeRule.newDoubleRule(STEP_SIZE),
            AttributeRule.newDoubleRule(KAPPA, true),
            new ElementRule(Parameter.class),
            new ElementRule(BayesianBridgeLikelihood.class),
            new ElementRule(GradientWrtParameterProvider.class, true),
            new ElementRule(SCORE_TERMS, new XMLSyntaxRule[]{
                    new ElementRule(Likelihood.class),
            }),
            new ElementRule(MASS_VARIANCE, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class),
            }),
    };
}
