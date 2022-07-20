/*
 * LaplaceApproximationParser.java
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

import dr.evomodel.treedatalikelihood.discrete.LaplaceApproximation;
import dr.evomodel.treedatalikelihood.discrete.MaximizerWrtParameter;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.util.Transform;
import dr.xml.*;

/**
 * @author Andy Magee
 */

public class LaplaceApproximationParser extends AbstractXMLObjectParser {

    private static final String NAME = "laplaceApproximation";
    private static final String DIAGONAL = "reportDiagonal";
    private static final String ML = "reportMarginalLikelihood";
    private static final String KL = "estimateKL";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        MaximizerWrtParameter maximizer = (MaximizerWrtParameter) xo.getChild(MaximizerWrtParameter.class);

        boolean diagonal = xo.getAttribute(DIAGONAL, true);
        boolean kl = xo.getAttribute(KL, true);
        boolean ml = xo.getAttribute(ML, true);

        LaplaceApproximation laplace = new LaplaceApproximation(maximizer, diagonal, kl, ml);

        laplace.approximate();

        return laplace;
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    @Override
    public String getParserDescription() {
        return null;
    }

    @Override
    public Class getReturnType() {
        return LaplaceApproximation.class;
    }

    @Override
    public String getParserName() {
        return NAME;
    }

    private static XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new ElementRule(MaximizerWrtParameter.class, false),
            AttributeRule.newBooleanRule(DIAGONAL, true, "Should we report only the diagonal of the covariance matrix?"),
            AttributeRule.newBooleanRule(KL, true, "Should the (variational-inference-style wrong-way) KL(aproximation||posterior) be estimated? Requires sampling from approximation."),
            AttributeRule.newBooleanRule(ML, true, "Should the (log) Marginal Likelihood (using the Laplace approximation) be reported?")
    };
}
