/*
 * ExpressionStatisticParser.java
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

import dr.inference.model.ExpressionStatistic;
import dr.inference.model.Statistic;
import dr.xml.*;

/**
 */
public class ExpressionStatisticParser extends AbstractXMLObjectParser {

    public static String EXPRESSION_STATISTIC = "expressionStatistic";
    public static String VARIABLES = "variables";
    public static String EXPRESSION = "expression";

    public String getParserName() {
        return EXPRESSION_STATISTIC;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        XMLObject cxo = xo.getChild(EXPRESSION);
        String expression = cxo.getStringChild(0);

        ExpressionStatistic expStatistic = new ExpressionStatistic(EXPRESSION_STATISTIC, expression);
        System.out.println("Expression: " + expression);
        System.out.println("  variables:");

        cxo = xo.getChild(VARIABLES);
        for (int i = 0; i < cxo.getChildCount(); i++) {
            Object child = cxo.getChild(i);
            if (child instanceof Statistic) {
                Statistic stat = (Statistic) child;

                expStatistic.addStatistic(stat);
                System.out.println("    " + stat.getStatisticName() + "=" + stat.getStatisticValue(0));
            } else {
                throw new XMLParseException("Unknown element found in " + getParserName() + " element:" + child);
            }
        }

        System.out.println("  value: " + expStatistic.getStatisticValue(0));


        return expStatistic;
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element returns a statistic that is the mean of the child statistics.";
    }

    public Class getReturnType() {
        return ExpressionStatistic.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new ElementRule(EXPRESSION,
                    new XMLSyntaxRule[]{new ElementRule(String.class)}),
            new ElementRule(VARIABLES,
                    new XMLSyntaxRule[]{new ElementRule(Statistic.class, 1, Integer.MAX_VALUE)})
    };
}
