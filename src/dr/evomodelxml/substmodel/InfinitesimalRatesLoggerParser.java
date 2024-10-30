/*
 * InfinitesimalRatesLoggerParser.java
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

package dr.evomodelxml.substmodel;

import dr.evomodel.substmodel.InfinitesimalRatesLogger;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.xml.*;

public class InfinitesimalRatesLoggerParser extends AbstractXMLObjectParser {

    private static final String NAME = "infinitesimalRatesLogger";
    private static final String DIAGONAL_ELEMENTS = "diagonalElements";
    private static final String LOGTRANSFORM = "logTransform";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        boolean diagonalElements = true;
        if (xo.hasAttribute(DIAGONAL_ELEMENTS)) {
            diagonalElements = xo.getAttribute(DIAGONAL_ELEMENTS, true);
        }

        boolean logTransform = false;
        if (xo.hasAttribute(LOGTRANSFORM)) {
            logTransform = xo.getAttribute(LOGTRANSFORM, true);
            if (diagonalElements && logTransform) {
                throw new XMLParseException("Cannot log-transform (negative) diagonal elements of the infinitesimal matrix");
            }
        }

        SubstitutionModel substitutionModel = (SubstitutionModel) xo.getChild(SubstitutionModel.class);
        return new InfinitesimalRatesLogger(substitutionModel, diagonalElements, logTransform);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        AttributeRule.newBooleanRule(DIAGONAL_ELEMENTS, true);
        AttributeRule.newBooleanRule(LOGTRANSFORM, true);
        return new XMLSyntaxRule[] {
                new ElementRule(SubstitutionModel.class),
        };
    }

    @Override
    public String getParserDescription() {
        return "Logger to report infinitesimal rates of a substitution model";
    }

    @Override
    public Class getReturnType() {
        return InfinitesimalRatesLogger.class;
    }

    @Override
    public String getParserName() {
        return NAME;
    }
}
