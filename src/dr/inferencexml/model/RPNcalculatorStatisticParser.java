/*
 * RPNcalculatorStatisticParser.java
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

import dr.inference.model.RPNcalculatorStatistic;
import dr.inference.model.Statistic;
import dr.xml.*;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 */
public class RPNcalculatorStatisticParser extends AbstractXMLObjectParser {

    public static String RPN_STATISTIC = "RPNcalculator";
    public static String VARIABLE = "variable";
    public static String EXPRESSION = "expression";
    public String getParserName() { return RPN_STATISTIC; }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        List<String> expressions =  new ArrayList<String>();
        List<String> expressionNames = new ArrayList<String>();

        Map<String, Statistic> variables = new HashMap<String, Statistic>();
        for(int i = 0; i < xo.getChildCount(); ++i) {
            final XMLObject child = (XMLObject)xo.getChild(i);
            if( child.getName().equals(EXPRESSION) ) {
                expressions.add(child.getStringChild(0));

                String name = child.hasAttribute(Statistic.NAME) ? child.getStringAttribute(Statistic.NAME) :
                        ("expression_" + expressionNames.size());

                expressionNames.add(name);

            } else if( child.getName().equals(VARIABLE) ) {
                Statistic s = (Statistic)child.getChild(Statistic.class);

                if( s.getDimension() != 1 )  {
                     throw new XMLParseException("Sorry, no support for multi-dimentional yet");
                }
                //assert s.getDimension() == 1;   // for now

                String name = child.hasAttribute(Statistic.NAME) ? child.getStringAttribute(Statistic.NAME) :
                        s.getDimensionName(0);
                variables.put(name, s);
            } else {
                 throw new XMLParseException("Unknown element found in " + getParserName() + " element:" + child);
            }
        }

        final String name = xo.hasAttribute(Statistic.NAME) ? xo.getStringAttribute(Statistic.NAME) : RPN_STATISTIC;
        final String[] e = expressions.toArray(new String[expressions.size()]);
        final String[] enames = expressionNames.toArray(new String[expressionNames.size()]);
        try {
            return new RPNcalculatorStatistic(name, e, enames, variables);
        } catch ( RuntimeException err) {
            throw new XMLParseException(err.getMessage());
        }        
    }

    public String getParserDescription() {
        return "This element returns a statistic evaluated from arbitrary expression.";
    }

    public Class getReturnType() { return RPNcalculatorStatistic.class; }

    public XMLSyntaxRule[] getSyntaxRules() { return rules; }

    private final XMLSyntaxRule[] rules = {
        new ElementRule(EXPRESSION,
            new XMLSyntaxRule[] { new ElementRule(String.class) }, 1, Integer.MAX_VALUE),
        new ElementRule(VARIABLE,
            new XMLSyntaxRule[] { new ElementRule(Statistic.class) } , 1, Integer.MAX_VALUE)
    };

}
