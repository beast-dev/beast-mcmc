/*
 * ValuesPoolSwapOperatorParser.java
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

package dr.inferencexml.operators;

import dr.inference.model.ValuesPool;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.ValuesPoolSwapOperator;
import dr.xml.*;

/**
 *
 */
public class ValuesPoolSwapOperatorParser extends AbstractXMLObjectParser {
    public static String VALUESPOOL_OPERATOR = "poolSwapOperator";

    public String getParserName() {
        return VALUESPOOL_OPERATOR;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        final ValuesPool parameter = (ValuesPool) xo.getChild(ValuesPool.class);
        final double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);
        final MCMCOperator op = new ValuesPoolSwapOperator(parameter) ;
        op.setWeight(weight);
        return op;
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return ".";
    }

    public Class getReturnType() {
        return ValuesPoolSwapOperator.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[] {
                AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
                new ElementRule(ValuesPool.class),
        };
    }

}
