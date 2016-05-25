/*
 * ValuesPoolParser.java
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
public class ValuesPoolParser extends dr.xml.AbstractXMLObjectParser {
    public static String VALUES_POOL = "valuesPool";
    public static String VALUES = "values";
    public static String SELECTOR = "selector";
    private static final String DEFAULT_VALUE = "default";

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        final Variable<Double> pool = (Variable<Double>)xo.getElementFirstChild(VALUES);
        final Variable<Double> selector = (Variable<Double>)xo.getElementFirstChild(SELECTOR);
        final double d = xo.getDoubleAttribute(DEFAULT_VALUE);

        if( pool.getSize() != selector.getSize() ) {
            throw new XMLParseException("variable Pool and selector must have equal length.");
        }
        return new ValuesPool(pool, selector, d);
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[] {
                AttributeRule.newDoubleRule(DEFAULT_VALUE),
                new ElementRule(VALUES, new XMLSyntaxRule[]{
                        new ElementRule(Variable.class,1,1) }),
                new ElementRule(SELECTOR, new XMLSyntaxRule[]{
                        new ElementRule(Variable.class,1,1) })
          };
      }


    public String getParserDescription() {
        return "";
    }

    public Class getReturnType() {
        return ValuesPool.class;
    }

    public String getParserName() {
        return VALUES_POOL;
    }
}
