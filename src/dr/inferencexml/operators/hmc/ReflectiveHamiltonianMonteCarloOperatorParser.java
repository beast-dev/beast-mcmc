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
import dr.inference.model.GeneralBoundsProvider;
import dr.inference.model.GraphicalParameterBound;
import dr.inference.model.Parameter;
import dr.inference.operators.AdaptationMode;
import dr.inference.operators.hmc.HamiltonianMonteCarloOperator;
import dr.inference.operators.hmc.MassPreconditionScheduler;
import dr.inference.operators.hmc.MassPreconditioner;
import dr.inference.operators.hmc.ReflectiveHamiltonianMonteCarloOperator;
import dr.util.Transform;
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
    private GeneralBoundsProvider bounds;

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        this.bounds = (GeneralBoundsProvider) xo.getChild(GeneralBoundsProvider.class);
        return super.parseXMLObject(xo);
    }

    @Override
    protected HamiltonianMonteCarloOperator factory(AdaptationMode adaptationMode, double weight, GradientWrtParameterProvider derivative,
                                                    Parameter parameter, Transform transform, Parameter mask,
                                                    HamiltonianMonteCarloOperator.Options runtimeOptions,
                                                    MassPreconditioner preconditioner, MassPreconditionScheduler.Type schedulerType) {

        return new ReflectiveHamiltonianMonteCarloOperator(adaptationMode, weight, derivative,
                parameter, transform, mask,
                runtimeOptions, preconditioner, bounds);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        XMLSyntaxRule[] extendedRules = new XMLSyntaxRule[rules.length + 1];
        extendedRules[0] = new ElementRule(GeneralBoundsProvider.class);
        for (int i = 0; i < rules.length; i++) {
            extendedRules[i + 1] = rules[i];
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
