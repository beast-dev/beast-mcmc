/*
 * FullyConjugateTreeTipsPotentialDerivativeParser.java
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

package dr.evomodelxml.continuous.hmc;

import dr.evomodel.continuous.FullyConjugateMultivariateTraitLikelihood;
import dr.evomodel.continuous.hmc.FullyConjugateTreeTipsPotentialDerivative;
import dr.xml.*;

/**
 * @author Max Tolkoff
 */
public class FullyConjugateTreeTipsPotentialDerivativeParser extends AbstractXMLObjectParser {
    public final static String FULLY_CONJUGATE_TREE_TIPS_POTENTIAL_DERIVATIVE = "fullyConjugateTreeTipsPotentialDerivative";
    public final static String FULLY_CONJUGATE_TREE_TIPS_POTENTIAL_DERIVATIVE2 = "traitGradientOnTree";


    @Override
    public String getParserName() {
        return FULLY_CONJUGATE_TREE_TIPS_POTENTIAL_DERIVATIVE;
    }

    @Override
    public String[] getParserNames() {
        return new String[] { FULLY_CONJUGATE_TREE_TIPS_POTENTIAL_DERIVATIVE, FULLY_CONJUGATE_TREE_TIPS_POTENTIAL_DERIVATIVE2 };
    }

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        FullyConjugateMultivariateTraitLikelihood treeLikelihood = (FullyConjugateMultivariateTraitLikelihood) xo.getChild(FullyConjugateMultivariateTraitLikelihood.class);

        return new FullyConjugateTreeTipsPotentialDerivative(treeLikelihood);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(FullyConjugateMultivariateTraitLikelihood.class),
    };

    @Override
    public String getParserDescription() {
        return null;
    }

    @Override
    public Class getReturnType() {
        return FullyConjugateTreeTipsPotentialDerivative.class;
    }
}
