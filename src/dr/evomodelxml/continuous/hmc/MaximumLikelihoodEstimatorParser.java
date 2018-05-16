/*
 * MaximumLikelihoodEstimatorParser.java
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

import dr.evomodel.treedatalikelihood.discrete.MaximumLikelihoodEstimator;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

/**
 * @author Marc A. Suchard
 * @author Xiang Ji
 */

public class MaximumLikelihoodEstimatorParser extends AbstractXMLObjectParser {

    private static final String NAME = "maximumLikelihoodEstimator";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        final GradientWrtParameterProvider gradientWrtParameterProvider = (GradientWrtParameterProvider) xo.getChild(GradientWrtParameterProvider.class);
        return new MaximumLikelihoodEstimator(gradientWrtParameterProvider);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[0];
    }

    @Override
    public String getParserDescription() {
        return null;
    }

    @Override
    public Class getReturnType() {
        return MaximumLikelihoodEstimator.class;
    }

    @Override
    public String getParserName() {
        return NAME;
    }
}
