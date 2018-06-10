/*
 * PrecisionGradientParser.java
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

import dr.evomodel.treedatalikelihood.continuous.WishartStatisticsWrapper;
import dr.evomodel.treedatalikelihood.hmc.CorrelationPrecisionGradient;
import dr.evomodel.treedatalikelihood.hmc.DiagonalPrecisionGradient;
import dr.evomodel.treedatalikelihood.hmc.PrecisionGradient;
import dr.inference.model.Likelihood;
import dr.inference.model.MatrixParameterInterface;
import dr.xml.*;

/**
 * @author Paul Bastide
 * @author Marc A. Suchard
 */

public class PrecisionGradientParser extends AbstractXMLObjectParser {

    private final static String PRECISION_GRADIENT = "precisionGradient";
    private final static String PARAMETER = "parameter";
    private final static String PRECISION_CORRELATION = "precisioncorrelation";
    private final static String PRECISION_DIAGONAL = "precisiondiagonal";
    private final static String PRECISION_BOTH = "precision";

    @Override
    public String getParserName() {
        return PRECISION_GRADIENT;
    }

    private int parseParameterMode(XMLObject xo) throws XMLParseException {
        // Choose which parameter(s) to update:
        // 0: full precision
        // 1: only precision correlation
        // 2: only precision diagonal
        int mode = 0;
        String parameterString = xo.getAttribute(PARAMETER, PRECISION_BOTH).toLowerCase();
        if (parameterString.compareTo(PRECISION_CORRELATION) == 0) {
            mode = 1;
        } else if (parameterString.compareTo(PRECISION_DIAGONAL) == 0) {
            mode = 2;
        }
        return mode;
    }

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        WishartStatisticsWrapper wishartStatistics = (WishartStatisticsWrapper)
                xo.getChild(WishartStatisticsWrapper.class);

        Likelihood likelihood = (Likelihood) xo.getChild(Likelihood.class);

        MatrixParameterInterface parameter = (MatrixParameterInterface) xo.getChild(MatrixParameterInterface.class);

        int parameterMode = parseParameterMode(xo);
        if (parameterMode == 0) {
            return new PrecisionGradient(wishartStatistics, likelihood, parameter);
        } else if (parameterMode == 1) {
            return new CorrelationPrecisionGradient(wishartStatistics, likelihood, parameter);
        } else {
            return new DiagonalPrecisionGradient(wishartStatistics, likelihood, parameter);
        }
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(WishartStatisticsWrapper.class),
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
