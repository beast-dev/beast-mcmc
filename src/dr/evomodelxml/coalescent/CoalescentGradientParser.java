/*
 * CoalescentGradientParser.java
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


package dr.evomodelxml.coalescent;

import dr.evolution.coalescent.CoalescentGradient;
import dr.evomodel.coalescent.CoalescentLikelihood;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Parameter;
import dr.inferencexml.operators.hmc.HamiltonianMonteCarloOperatorParser;
import dr.xml.*;

import static dr.evolution.coalescent.CoalescentGradient.Wrt;

/**
 * @author Xiang Ji
 * @author Marc A. Suchard
 */
public class CoalescentGradientParser extends AbstractXMLObjectParser {

    private static final String NAME = "coalescentGradient";

    private static final String TOLERANCE = HamiltonianMonteCarloOperatorParser.GRADIENT_CHECK_TOLERANCE;
    private static final String WRT = "wrt";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        double tolerance = xo.getAttribute(TOLERANCE, 1E-1);
        CoalescentLikelihood likelihood = (CoalescentLikelihood) xo.getChild(CoalescentLikelihood.class);
        TreeModel tree = (TreeModel) xo.getChild(TreeModel.class);

        if (xo.hasChildNamed(WRT)) {
            Parameter wrtParameter = (Parameter) xo.getElementFirstChild(WRT);
            return new CoalescentGradient(likelihood, tree, wrtParameter, Wrt.PARAMETER, tolerance);
        } else {
            return new CoalescentGradient(likelihood, tree, null, Wrt.NODE_HEIGHTS, tolerance);
        }
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(CoalescentLikelihood.class),
            new ElementRule(TreeModel.class),
            new ElementRule(WRT, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class),
            }, true),
            AttributeRule.newDoubleRule(TOLERANCE, true),
    };


    @Override
    public String getParserDescription() {
        return null;
    }

    @Override
    public Class getReturnType() {
        return CoalescentGradient.class;
    }

    @Override
    public String getParserName() {
        return NAME;
    }
}
