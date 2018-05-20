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

import dr.evomodel.treedatalikelihood.hmc.PrecisionGradient;
import dr.inference.model.Likelihood;
import dr.inference.model.MatrixParameterInterface;
import dr.inference.model.Parameter;
import dr.math.interfaces.ConjugateWishartStatisticsProvider;
import dr.xml.*;

/**
 * @author Paul Bastide
 * @author Marc A. Suchard
 */

public class PrecisionGradientParser extends AbstractXMLObjectParser {

    private final static String PRECISION_GRADIENT = "precisionGradient";

    @Override
    public String getParserName() {
        return PRECISION_GRADIENT;
    }


    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        ConjugateWishartStatisticsProvider wishartStatistics = (ConjugateWishartStatisticsProvider)
                xo.getChild(ConjugateWishartStatisticsProvider.class);

        Likelihood likelihood = (Likelihood) xo.getChild(Likelihood.class);

        // TODO Get information about precision parametrization
        Parameter parameter = (MatrixParameterInterface) xo.getChild(MatrixParameterInterface.class);

        return new PrecisionGradient(wishartStatistics, likelihood, parameter);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(ConjugateWishartStatisticsProvider.class),
            new ElementRule(Likelihood.class),
            new ElementRule(MatrixParameterInterface.class),
    };

    @Override
    public String getParserDescription() {
        return null;
    }

    @Override
    public Class getReturnType() {
        return PrecisionGradient.class;
    }
}
