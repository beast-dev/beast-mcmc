/*
 * DiagonalPriorPreconditionerParser.java
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

package dr.inferencexml.hmc;

import dr.inference.hmc.DiagonalPriorPreconditioner;
import dr.inference.model.Parameter;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

/**
 * Parses a fixed vector of prior SDs for HMC prior-diagonal preconditioning.
 *
 * @author Filippo Monti
 */
public class DiagonalPriorPreconditionerParser extends AbstractXMLObjectParser {

    public static final String DIAGONAL_PRIOR_PRECONDITIONER = "diagonalPriorPreconditioner";
    private static final String STANDARD_DEVIATION = "standardDeviation";

    @Override
    public Object parseXMLObject(final XMLObject xo) throws XMLParseException {
        final Parameter sd = (Parameter) xo.getElementFirstChild(STANDARD_DEVIATION);
        final String id = xo.hasId() ? xo.getId() : DIAGONAL_PRIOR_PRECONDITIONER;
        return new DiagonalPriorPreconditioner(id, sd.getParameterValues());
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(STANDARD_DEVIATION, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }),
    };

    @Override
    public String getParserDescription() {
        return "A fixed vector of prior standard deviations for HMC prior-diagonal preconditioning.";
    }

    @Override
    public Class getReturnType() {
        return DiagonalPriorPreconditioner.class;
    }

    @Override
    public String getParserName() {
        return DIAGONAL_PRIOR_PRECONDITIONER;
    }
}
