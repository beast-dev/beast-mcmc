/*
 * ReflectiveHamiltonianMonteCarloOperatorParser.java
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

import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.GeneralBoundsProvider;
import dr.inference.model.GraphicalParameterBound;
import dr.inference.model.Parameter;
import dr.inference.operators.AdaptationMode;
import dr.inference.operators.hmc.HamiltonianMonteCarloOperator;
import dr.inference.operators.hmc.MassPreconditionScheduler;
import dr.inference.operators.hmc.MassPreconditioner;
import dr.inference.operators.hmc.ReflectiveHamiltonianMonteCarloOperator;
import dr.util.Transform;
import dr.xml.AttributeRule;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

/**
 * @author Xiang Ji
 * @author Marc A. Suchard
 */
public class ReflectiveHamiltonianMonteCarloOperatorParser extends HamiltonianMonteCarloOperatorParser {

    public final static String OPERATOR_NAME = "reflectiveHamiltonianMonteCarloOperator";
    public final static String REFLECTION_RULE = "reflectionRule";
    public final static String MAX_EVENTS_PER_STEP = "maxEventsPerStep";
    private GeneralBoundsProvider bounds;
    private ReflectiveHamiltonianMonteCarloOperator.ReflectionRule reflectionRule;
    private long maxEventsPerStep;

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        if (xo.hasChildNamed(MASK)) {
            throw new XMLParseException(OPERATOR_NAME + " does not support a mask; remove the <" + MASK + "> element");
        }
        this.bounds = (GeneralBoundsProvider) xo.getChild(GeneralBoundsProvider.class);
        try {
            this.reflectionRule = ReflectiveHamiltonianMonteCarloOperator.ReflectionRule.parse(xo.getAttribute(
                    REFLECTION_RULE, ReflectiveHamiltonianMonteCarloOperator.DEFAULT_REFLECTION_RULE.getName()));
        } catch (IllegalArgumentException e) {
            throw new XMLParseException(e.getMessage());
        }
        final double maxEvents = xo.getAttribute(MAX_EVENTS_PER_STEP,
                (double) ReflectiveHamiltonianMonteCarloOperator.DEFAULT_MAX_EVENTS_PER_STEP);
        if (!(maxEvents >= 1) || maxEvents > Long.MAX_VALUE || maxEvents != Math.floor(maxEvents)) {
            throw new XMLParseException(MAX_EVENTS_PER_STEP + " must be a whole number of at least 1");
        }
        this.maxEventsPerStep = (long) maxEvents;
        if (bounds instanceof GraphicalParameterBound) {
            // the operator's parameter, found as HamiltonianMonteCarloOperatorParser.parseXMLObject finds it
            Parameter parameter = (Parameter) xo.getChild(Parameter.class);
            if (parameter == null) {
                parameter = ((GradientWrtParameterProvider) xo.getChild(GradientWrtParameterProvider.class))
                        .getParameter();
            }
            try {
                ((GraphicalParameterBound) bounds).checkParameterMapping();
                ((GraphicalParameterBound) bounds).checkParameterMapping(parameter);
            } catch (IllegalArgumentException e) {
                throw new XMLParseException(e.getMessage());
            }
        }
        return super.parseXMLObject(xo);
    }

    @Override
    protected HamiltonianMonteCarloOperator factory(AdaptationMode adaptationMode, double weight, GradientWrtParameterProvider derivative,
                                                    Parameter parameter, Transform transform, Parameter mask,
                                                    HamiltonianMonteCarloOperator.Options runtimeOptions,
                                                    MassPreconditioner preconditioner, MassPreconditionScheduler.Type schedulerType) {

        return new ReflectiveHamiltonianMonteCarloOperator(adaptationMode, weight, derivative,
                parameter, transform, mask,
                runtimeOptions, preconditioner, bounds, reflectionRule, maxEventsPerStep);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        XMLSyntaxRule[] extendedRules = new XMLSyntaxRule[rules.length + 3];
        extendedRules[0] = new ElementRule(GeneralBoundsProvider.class);
        extendedRules[1] = AttributeRule.newStringRule(REFLECTION_RULE, true);
        extendedRules[2] = AttributeRule.newDoubleRule(MAX_EVENTS_PER_STEP, true);
        for (int i = 0; i < rules.length; i++) {
            extendedRules[i + 3] = rules[i];
        }
        return extendedRules;
    }


    @Override
    public String getParserDescription() {
        return "Returns a Reflective Hamiltonian Monte Carlo transition kernel";
    }

    @Override
    public Class getReturnType() {
        return ReflectiveHamiltonianMonteCarloOperator.class;
    }

    @Override
    public String getParserName() {
        return OPERATOR_NAME;
    }
}
