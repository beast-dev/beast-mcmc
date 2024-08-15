/*
 * BoundedSpaceCorrelationParser.java
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

package dr.inferencexml.model;

import dr.inference.model.BoundedSpace;
import dr.xml.*;

public class BoundedSpaceCorrelationParser extends AbstractXMLObjectParser {
    private static final String CORRELATION_BOUNDS = "correlationBounds";
    private static final String DIMENSION = "dimension";


    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        int dim = xo.getIntegerAttribute(DIMENSION);
        return new BoundedSpace.Correlation(dim);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                AttributeRule.newIntegerRule(DIMENSION)
        };
    }

    @Override
    public String getParserDescription() {
        return "Indicates whether a parameter is a valid correlation matrix or not";
    }

    @Override
    public Class getReturnType() {
        return BoundedSpace.Correlation.class;
    }

    @Override
    public String getParserName() {
        return CORRELATION_BOUNDS;
    }
}
