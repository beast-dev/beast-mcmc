/*
 * BitFlipOperator.java
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

package dr.inference.operators;

import dr.inference.model.Parameter;
import dr.math.MathUtils;
import dr.xml.*;

/**
 * A generic operator that flips bits.
 *
 * @author Alexei Drummond
 *
 * @version $Id: DeltaExchangeOperator.java,v 1.18 2005/06/14 10:40:34 rambaut Exp $
 */
public class BitFlipOperator extends SimpleMCMCOperator {
	
	public static final String BIT_FLIP_OPERATOR = "bitFlipOperator";

	public BitFlipOperator(Parameter parameter, int weight) {
		this.parameter = parameter;
		this.weight = weight;
	}
	
	/** @return the parameter this operator acts on. */
	public Parameter getParameter() { return parameter; }

	/** change the parameter and return the hastings ratio. 
	 * performs a delta exchange operation between two scalars in the vector
	 * and return the hastings ratio.
	 */
	public final double doOperation() throws OperatorFailedException {


        int n = parameter.getDimension();
        double sum = 0.0;
        double logq = 0.0;

        for (int i = 0; i < n; i++) {
            sum += parameter.getParameterValue(i);
        }


        int pos = MathUtils.nextInt(parameter.getDimension());

        int value = (int)parameter.getParameterValue(pos);

        if (value == 0) {
            parameter.setParameterValue(pos, 1.0);

            logq = -Math.log((n-sum) / (sum+1));

        } else if (value == 1) {
            parameter.setParameterValue(pos, 0.0);

            logq = -Math.log(sum / (n-sum-1) );

        } else {
            throw new RuntimeException("expected 1 or 0");
        }

        // hastings ratio is designed to make move symmetric on sum of 1's
        return logq;
	}

	// Interface MCMCOperator
	public final String getOperatorName() { return parameter.getParameterName(); }

	public int getWeight() { return weight; }

	public void setWeight(int w) { weight = w; }

	public final String getPerformanceSuggestion() {

		return "no performance suggestion";
	}
	
	public String toString() { return getOperatorName(); }
	
	public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {
		
		public String getParserName() { return BIT_FLIP_OPERATOR; }
		
		public Object parseXMLObject(XMLObject xo) throws XMLParseException {

			int weight = xo.getIntegerAttribute(WEIGHT);

			Parameter parameter = (Parameter)xo.getChild(Parameter.class);	
			
			return new BitFlipOperator(parameter, weight);
		}
		
		//************************************************************************
		// AbstractXMLObjectParser implementation
		//************************************************************************
		
		public String getParserDescription() {
			return "This element returns a bit-flip operator on a given parameter.";
		}
		
		public Class getReturnType() { return MCMCOperator.class; }
		
		public XMLSyntaxRule[] getSyntaxRules() { return rules; }
		
		private XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
			AttributeRule.newIntegerRule(WEIGHT),
			new ElementRule(Parameter.class)
		};
	
	};
	// Private instance variables
	
	private Parameter parameter = null;
	private int weight = 1;
}
