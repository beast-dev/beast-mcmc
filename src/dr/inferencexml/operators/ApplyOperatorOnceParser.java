/*
 * ApplyOperatorOnceParser.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.inferencexml.operators;

import dr.inference.operators.SimpleMCMCOperator;
import dr.xml.*;

/**
 * @author Marc Suchard
 * @author Mandev Gill
 */
public class ApplyOperatorOnceParser extends AbstractXMLObjectParser {

    public static final String APPLY = "applyOperatorOnce";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        for (int i = 0; i < xo.getChildCount(); ++i) {
            SimpleMCMCOperator operator = (SimpleMCMCOperator) xo.getChild(i);
            operator.doOperation();
        }
        return null;
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                new ElementRule(SimpleMCMCOperator.class, 1, Integer.MAX_VALUE),
        };
    }

    @Override
    public String getParserDescription() {
        return "Applies a series of operators once before continuing";
    }

    @Override
    public Class getReturnType() {
        return ApplyOperatorOnceParser.class;
    }

    public String getParserName() {
        return APPLY;
    }
}
