/*
 * ParameterChooserParser.java
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

package dr.inference.model;

import dr.xml.*;

/**
 * @author Joseph Heled
 *         Date: 4/09/2009
 */
public class ParameterChooserParser extends dr.xml.AbstractXMLObjectParser {
    public static String VARIABLE_SELECTOR = "variableSelector";
    public static String INDEX = "index";

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        String name = xo.getId();
        final ValuesPool pool = (ValuesPool) xo.getChild(ValuesPool.class);
        final int[] which = xo.getIntegerArrayAttribute(INDEX);
        for( int w : which ) {
            if( ! (0 <= w && w < pool.length()) ) {
                throw new XMLParseException("index " + w + " out of range");
            }
        }
        return new ParameterChooser(name, pool, which);
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[] {
                new ElementRule(ValuesPool.class,1,1),
                AttributeRule.newIntegerArrayRule(INDEX, false),
        };
    }

    public String getParserDescription() {
        return "";
    }

    public Class getReturnType() {
        return ParameterChooser.class;
    }

    public String getParserName() {
        return VARIABLE_SELECTOR;
    }
}
