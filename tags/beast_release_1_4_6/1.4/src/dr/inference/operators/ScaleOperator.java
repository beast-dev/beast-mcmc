/*
 * ScaleOperator.java
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
 * A generic scale operator for use with a multi-dimensional parameters.
 *
 * @author Alexei Drummond
 * @author Andrew Rambaut
 *
 * @version $Id: ScaleOperator.java,v 1.20 2005/06/14 10:40:34 rambaut Exp $
 */
public class ScaleOperator extends SimpleMCMCOperator implements CoercableMCMCOperator {
	
	public static final String SCALE_OPERATOR = "scaleOperator";
	public static final String SCALE_ALL = "scaleAll";
	public static final String SCALE_FACTOR = "scaleFactor";

	public ScaleOperator(Parameter parameter, boolean scaleAll, double scale, int weight, int mode) {
		this.parameter = parameter;
		this.scaleAll = scaleAll;
		this.scaleFactor = scale;
		this.weight = weight;
		this.mode = mode;
	}
	
	/** @return the parameter this operator acts on. */
	public Parameter getParameter() { return parameter; }

	/**
	 * change the parameter and return the hastings ratio.
	 */
	public final double doOperation() throws OperatorFailedException {
		
		double scale = (scaleFactor + (MathUtils.nextDouble() * ((1.0/scaleFactor) - scaleFactor)));

		double logq;
		
		if (scaleAll) {
			// update all dimensions

			logq = (parameter.getDimension() - 2) * Math.log(scale);

			for (int i = 0; i < parameter.getDimension(); i++) {
				parameter.setParameterValue(i, parameter.getParameterValue(i) * scale);
			}
			for (int i =0; i < parameter.getDimension(); i++) {
				if (parameter.getParameterValue(i) < parameter.getBounds().getLowerLimit(i) || 
					parameter.getParameterValue(i) > parameter.getBounds().getUpperLimit(i)) {
					throw new OperatorFailedException("proposed value outside boundaries");
				}
			}
		} else {
    		logq = - Math.log(scale);

			int index = MathUtils.nextInt(parameter.getDimension());
			
			double oldValue = parameter.getParameterValue(index);
            double newValue = scale * oldValue;
				
			if (newValue < parameter.getBounds().getLowerLimit(index) || 
				newValue > parameter.getBounds().getUpperLimit(index)) {
				throw new OperatorFailedException("proposed value outside boundaries");
			}
			
			parameter.setParameterValue(index, newValue);

            // provides a hook for subclasses
            cleanupOperation(newValue, oldValue);
		}

		return logq;
	}

    /**
     * This method should be overridden by operators that need to do something just before the return of doOperation.
     * @param newValue the proposed parameter value
     * @param oldValue the old parameter value
     */
    void cleanupOperation(double newValue, double oldValue) {
        // DO NOTHING
    }

	//MCMCOperator INTERFACE
	public final String getOperatorName() { return parameter.getParameterName(); }

	public double getCoercableParameter() {
		return Math.log(1.0/scaleFactor - 1.0);
	}
	
	public void setCoercableParameter(double value) {
		scaleFactor = 1.0/(Math.exp(value) + 1.0);
	}
	
	public double getRawParameter() {
		return scaleFactor;
	}
	
	public int getMode() { 
		return mode; 
	}
	
	public double getScaleFactor() { return scaleFactor; }

	public double getTargetAcceptanceProbability() { return 0.234; }

	public double getMinimumAcceptanceLevel() { return 0.1;}
	public double getMaximumAcceptanceLevel() { return 0.4;}
	public double getMinimumGoodAcceptanceLevel() { return 0.20; }
	public double getMaximumGoodAcceptanceLevel() { return 0.30; }
	
	public int getWeight() { return weight; }
	public void setWeight(int w) { weight = w; }
	
	public final String getPerformanceSuggestion() {

		double prob = MCMCOperator.Utils.getAcceptanceProbability(this);
		double targetProb = getTargetAcceptanceProbability();
		dr.util.NumberFormatter formatter = new dr.util.NumberFormatter(5);
		double sf = OperatorUtils.optimizeScaleFactor(scaleFactor, prob, targetProb);
		if (prob < getMinimumGoodAcceptanceLevel()) {
			return "Try setting scaleFactor to about " + formatter.format(sf);
		} else if (prob > getMaximumGoodAcceptanceLevel()) {
			return "Try setting scaleFactor to about " + formatter.format(sf);
		} else return "";
	}
	
	public static dr.xml.XMLObjectParser PARSER = new dr.xml.AbstractXMLObjectParser() {
		
		public String getParserName() { return SCALE_OPERATOR; }
		
		public Object parseXMLObject(XMLObject xo) throws XMLParseException {
		
			boolean scaleAll = false;
			int mode = CoercableMCMCOperator.DEFAULT;

			
			if (xo.hasAttribute(SCALE_ALL)) {
				scaleAll = xo.getBooleanAttribute(SCALE_ALL);
			}
			
			if (xo.hasAttribute(AUTO_OPTIMIZE)) {
				if (xo.getBooleanAttribute(AUTO_OPTIMIZE)) {
					mode = CoercableMCMCOperator.COERCION_ON;
				} else {
					mode = CoercableMCMCOperator.COERCION_OFF;
				}
			}
			
			int weight = xo.getIntegerAttribute(WEIGHT);
			double scaleFactor = xo.getDoubleAttribute(SCALE_FACTOR);
			
			if (scaleFactor <= 0.0 || scaleFactor >= 1.0) {
				throw new XMLParseException("scaleFactor must be between 0.0 and 1.0");
			}
			
			Parameter parameter = (Parameter)xo.getChild(Parameter.class);	
			
			return new ScaleOperator(parameter, scaleAll, scaleFactor, weight, mode);
		}
		
		//************************************************************************
		// AbstractXMLObjectParser implementation
		//************************************************************************
		
		public String getParserDescription() {
			return "This element returns a scale operator on a given parameter.";
		}
		
		public Class getReturnType() { return MCMCOperator.class; }
		
		public XMLSyntaxRule[] getSyntaxRules() { return rules; }
		
		private XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
			AttributeRule.newDoubleRule(SCALE_FACTOR),
			AttributeRule.newBooleanRule(SCALE_ALL, true),
			AttributeRule.newIntegerRule(WEIGHT),
			AttributeRule.newBooleanRule(AUTO_OPTIMIZE, true),
			new ElementRule(Parameter.class)
		};
	
	};


	public String toString() { 
		return "scaleOperator(" + parameter.getParameterName() + " [" + scaleFactor + ", " + (1.0/scaleFactor) + "]"; 
	}
	
	//PRIVATE STUFF
	
	private Parameter parameter = null;
	private boolean scaleAll = false;
	private double scaleFactor = 0.5;
	private int mode = CoercableMCMCOperator.DEFAULT;
	private int weight = 1;
	
}
