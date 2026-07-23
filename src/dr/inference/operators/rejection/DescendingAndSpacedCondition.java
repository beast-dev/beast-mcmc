/*
 * DescendingAndSpacedCondition.java
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

package dr.inference.operators.rejection;

import dr.xml.*;

public class DescendingAndSpacedCondition implements AcceptCondition {
    private final double spacing;

    DescendingAndSpacedCondition(double spacing) {
        this.spacing = spacing;
    }

    @Override
    public boolean satisfiesCondition(double[] values) {
        for (int i = 1; i < values.length; i++) {
            if (spacing * Math.abs(values[i - 1]) < Math.abs(values[i])) {
                return false;
            }
        }
        return true;
    }

    private static final String DESCENDING_AND_SPACED = "descendingAndSpaced";
    private static final String SPACING = "spacing";

    public static AbstractXMLObjectParser PARSER = new AbstractXMLObjectParser() {
        @Override
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            double spacing = xo.getDoubleAttribute(SPACING);
            if (spacing < 0.0 || spacing > 1.0) {
                throw new XMLParseException("Attribute '" + SPACING + "' must be between 0 and 1.");
            }

            return new DescendingAndSpacedCondition(spacing);
        }

        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            return new XMLSyntaxRule[]{
                    AttributeRule.newDoubleRule(SPACING)
            };
        }

        @Override
        public String getParserDescription() {
            return "Condition requiring parameter to have descending absolute values with some minimum spacing.";
        }

        @Override
        public Class getReturnType() {
            return DescendingAndSpacedCondition.class;
        }

        @Override
        public String getParserName() {
            return DESCENDING_AND_SPACED;
        }
    };
}

