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
// **REMOVED**import org.nfunk.jep.JEP;

import java.util.Vector;

/**
 * @version $Id: ExpressionStatistic.java,v 1.4 2005/05/24 20:26:00 rambaut Exp $
 *
 * @author Alexei Drummond
 */
public class ExpressionStatistic extends Statistic.Abstract {
	
	public static String EXPRESSION_STATISTIC = "expressionStatistic";
    public static String VARIABLES = "variables";
    public static String EXPRESSION = "expression";

    String expression = "";

	public ExpressionStatistic(String name, String expression) {
		super(name);
        this.expression = expression;
	}
	
	public void addStatistic(Statistic statistic) {
        if (statistic.getDimension() != 1) {
            throw new IllegalArgumentException("Can only have statistics of dimension 1");
        }

		statistics.add(statistic);
	}
	
	public int getDimension() { return 1; }

	/** @return the value of the expression */
	public double getStatisticValue(int dim) {

		System.err.println("Error in parsing expression " + expression + " : JEP expression parser not included with this version");
        return 0;
	}
	
	
	public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {
		
		public String getParserName() { return EXPRESSION_STATISTIC; }
		
		public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            XMLObject cxo = (XMLObject)xo.getChild(EXPRESSION);
            String expression = cxo.getStringChild(0);

            ExpressionStatistic expStatistic = new ExpressionStatistic(EXPRESSION_STATISTIC, expression);
			System.out.println("Expression: " + expression);
            System.out.println("  variables:" );

			cxo = (XMLObject)xo.getChild(VARIABLES);
            for (int i =0; i < cxo.getChildCount(); i++) {
				Object child = cxo.getChild(i);
				if (child instanceof Statistic) {
                    Statistic stat = (Statistic)child;

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
		
		public Class getReturnType() { return ExpressionStatistic.class; }
		
		public XMLSyntaxRule[] getSyntaxRules() { return rules; }
		
		private XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
			new ElementRule(EXPRESSION,
				new XMLSyntaxRule[] { new ElementRule(String.class) }),
            new ElementRule(VARIABLES,
				new XMLSyntaxRule[] { new ElementRule(Statistic.class, 1, Integer.MAX_VALUE ) })
		};		
	};
	

	// ****************************************************************
	// Private and protected stuff
	// ****************************************************************
	
	private Vector statistics = new Vector();
}
