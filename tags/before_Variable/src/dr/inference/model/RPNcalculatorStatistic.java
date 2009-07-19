/*
 * ExpressionStatistic.java
 *
 * Copyright (C) 2002-2006 Alexei Drummond and Andrew Rambaut
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

import java.util.*;

/**
 *  A statistic based on evaluating simple expressions.
 *
 * The expressions are in RPN, so no parsing issues. whitspace separated. Variables (other statistics),
 * constants and operations. Currently just the basic four, but easy to extend.
 *
 * @author Joseph Heled
 */
public class RPNcalculatorStatistic extends Statistic.Abstract {

	public static String RPN_STATISTIC = "RPNcalculator";
    public static String VARIABLE = "variable";
    public static String EXPRESSION = "expression";

    private RPNexpressionCalculator[] expressions;
    private String[] names;
    private Map<String, Statistic> variables;

    RPNexpressionCalculator.GetVariable vars = new RPNexpressionCalculator.GetVariable() {
        public double get(String name) {
            return variables.get(name).getStatisticValue(0);
        }
    };

    public RPNcalculatorStatistic(String name, String[] expressions, String[] names,
                                  Map<String, Statistic> variables) {
		super(name);

        this.expressions = new RPNexpressionCalculator[expressions.length];
        for(int i = 0; i < expressions.length; ++ i) {
           this.expressions[i] = new RPNexpressionCalculator(expressions[i]);
        }

        this.names = names;
        this.variables = variables;
    }

	public int getDimension() {
        return expressions.length;
    }

    public String getDimensionName(int dim) {
        return names[dim];
    }

    /** @return the value of the expression */
	public double getStatisticValue(int dim) {
        return expressions[dim].evaluate(vars);
	}

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

		public String getParserName() { return RPN_STATISTIC; }

		public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            List<String> expressions =  new ArrayList<String>();
            List<String> expressionNames = new ArrayList<String>();

            Map<String, Statistic> variables = new HashMap<String, Statistic>();
            for(int i = 0; i < xo.getChildCount(); ++i) {
                final XMLObject child = (XMLObject)xo.getChild(i);
                if( child.getName().equals(EXPRESSION) ) {
                    expressions.add(child.getStringChild(0));
                    
                    String name = child.hasAttribute(NAME) ? child.getStringAttribute(NAME) :
                            ("expression_" + expressionNames.size());

                    expressionNames.add(name);

                } else if( child.getName().equals(VARIABLE) ) {
                    Statistic s = (Statistic)child.getChild(Statistic.class);

                    if( s.getDimension() != 1 )  {
                         throw new XMLParseException("Sorry, no support for multi-dimentional yet");
                    }
                    //assert s.getDimension() == 1;   // for now

                    String name = child.hasAttribute(NAME) ? child.getStringAttribute(NAME) :
                            s.getDimensionName(0);
                    variables.put(name, s);
                } else {
                     throw new XMLParseException("Unknown element found in " + getParserName() + " element:" + child);
                }
            }

            final String name = xo.hasAttribute(NAME) ? xo.getStringAttribute(NAME) : RPN_STATISTIC;
            final String[] e = expressions.toArray(new String[expressions.size()]);
            final String[] enames = expressionNames.toArray(new String[expressionNames.size()]);
            return new RPNcalculatorStatistic(name, e, enames, variables);
		}

		public String getParserDescription() {
			return "This element returns a statistic evaluated from arbitrary expression.";
		}

		public Class getReturnType() { return RPNcalculatorStatistic.class; }

		public XMLSyntaxRule[] getSyntaxRules() { return rules; }

		private XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
			new ElementRule(EXPRESSION,
				new XMLSyntaxRule[] { new ElementRule(String.class) }, 1, Integer.MAX_VALUE),
            new ElementRule(VARIABLE,
				new XMLSyntaxRule[] { new ElementRule(Statistic.class) } , 1, Integer.MAX_VALUE)
		};
	};
}