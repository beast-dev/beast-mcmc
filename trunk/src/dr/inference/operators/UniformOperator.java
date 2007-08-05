/*
 * UniformOperator.java
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
 * A generic uniform sampler/operator for use with a multi-dimensional parameter.
 *
 * @author Alexei Drummond
 * @author Andrew Rambaut
 *
 * @version $Id: UniformOperator.java,v 1.16 2005/06/14 10:40:34 rambaut Exp $
 */
public class UniformOperator extends SimpleMCMCOperator {
	
	public UniformOperator(Parameter parameter, int weight) {
		this.parameter = parameter;
		this.weight = weight;
	}
	
	/** @return the parameter this operator acts on. */
	public Parameter getParameter() { return parameter; }

	/**
	 * change the parameter and return the hastings ratio.
	 */
	public final double doOperation() {
		
		int index = MathUtils.nextInt(parameter.getDimension());
		double lower = parameter.getBounds().getLowerLimit(index);
		double upper = parameter.getBounds().getUpperLimit(index);
		double newValue = (MathUtils.nextDouble() * (upper - lower)) + lower;
		
		parameter.setParameterValue(index, newValue);
		
//		System.out.println(newValue + "[" + lower + "," + upper + "]");
		return 0.0;
	}

	//MCMCOperator INTERFACE
	public final String getOperatorName() {
        return "uniform(" + parameter.getParameterName() + ")"; 
    }

	public final void optimize(double targetProb) {	
	
		throw new RuntimeException("This operator cannot be optimized!");
	}

	public boolean isOptimizing() { return false; }
	public void setOptimizing(boolean opt) { throw new RuntimeException("This operator cannot be optimized!"); }

	public double getMinimumAcceptanceLevel() { return 0.1;}
	public double getMaximumAcceptanceLevel() { return 0.4;}
	public double getMinimumGoodAcceptanceLevel() { return 0.20; }
	public double getMaximumGoodAcceptanceLevel() { return 0.30; }

	public int getWeight() { return weight; }

	public void setWeight(int w) { weight = w; }
	
	public String getPerformanceSuggestion() { 
		if (MCMCOperator.Utils.getAcceptanceProbability(this) < getMinimumAcceptanceLevel()) {
			return "";
		} else if (MCMCOperator.Utils.getAcceptanceProbability(this) > getMaximumAcceptanceLevel()){
			return "";
		} else {
			return "";
		}
	}
	
	public static dr.xml.XMLObjectParser PARSER = new dr.xml.AbstractXMLObjectParser() {
		
		public String getParserName() { return "uniformOperator"; }
		
		public Object parseXMLObject(XMLObject xo) throws XMLParseException {
		
			int weight = xo.getIntegerAttribute(WEIGHT);
			Parameter parameter = (Parameter)xo.getChild(Parameter.class);	
			
			return new UniformOperator(parameter,weight);
		}
		
		//************************************************************************
		// AbstractXMLObjectParser implementation
		//************************************************************************
		
		public String getParserDescription() {
			return "An operator that picks new parameter values uniformly at random."; }
		
		public Class getReturnType() { return UniformOperator.class; }
		
		
		public XMLSyntaxRule[] getSyntaxRules() { return rules; }
		
		private XMLSyntaxRule[] rules = new XMLSyntaxRule[] {  
			AttributeRule.newIntegerRule(WEIGHT),
			new ElementRule(Parameter.class)
		};
	};

	public String toString() { return "uniformOperator(" + parameter.getParameterName() + ")"; }
	
	//PRIVATE STUFF
	
	private Parameter parameter = null;
	private int weight = 1;
}
