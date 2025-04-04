/*
 * TransformedParameterOperatorParser.java
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

import dr.inference.model.BoundedSpace;
import dr.inference.operators.SimpleMCMCOperator;
import dr.inference.operators.TransformedParameterOperator;
import dr.xml.*;

import static dr.inference.operators.TransformedParameterOperator.TRANSFORMED_OPERATOR;

public class TransformedParameterOperatorParser extends AbstractXMLObjectParser {


    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        SimpleMCMCOperator operator = (SimpleMCMCOperator) xo.getChild(SimpleMCMCOperator.class);
        BoundedSpace bounds = (BoundedSpace) xo.getChild(BoundedSpace.class);
        return new TransformedParameterOperator(operator, bounds);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                new ElementRule(SimpleMCMCOperator.class),
                new ElementRule(BoundedSpace.class, true)
        };
    }

    @Override
    public String getParserDescription() {
        return "operator that corrects the hastings ratio with appropriate Jacobian term due to parameter transform";
    }

    @Override
    public Class getReturnType() {
        return TransformedParameterOperator.class;
    }

    @Override
    public String getParserName() {
        return TRANSFORMED_OPERATOR;
    }
}
