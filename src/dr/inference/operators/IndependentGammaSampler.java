/*
 * IndependentGammaSampler.java
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

package dr.inference.operators;

import dr.inference.model.Bounds;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.math.distributions.GammaDistribution;
import dr.xml.AttributeRule;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

/**
 * An independent gamma sampler to propose new (independent) values from a provided gamma prior (distribution).
 *
 * @author Guy Baele
 * 
 */
public class IndependentGammaSampler extends SimpleMCMCOperator {

	public static final String OPERATOR_NAME = "independentGammaSampler";
    public static final String SHAPE = "shape";
    public static final String SCALE = "scale";
	
	private Variable<Double> variable = null;
	private GammaDistribution gamma = null;
	private boolean updateAllIndependently = true;
	
	public IndependentGammaSampler(Variable variable, GammaDistribution gamma) {
		
		this(variable, gamma, 1.0);
		
	}
	
	public IndependentGammaSampler(Variable variable, GammaDistribution gamma, double weight) {
		
		this(variable, gamma, weight, true);
		
	}
	
	public IndependentGammaSampler(Variable variable, GammaDistribution gamma, double weight, boolean updateAllIndependently) {
		
		this.variable = variable;
		this.gamma = gamma;
		this.updateAllIndependently = updateAllIndependently;
		setWeight(weight);
		
	}
	
	public String getPerformanceSuggestion() {
		return "";
	}

	public String getOperatorName() {
		return "independentGamma(" + variable.getVariableName() + ")";
	}

    /**
	 * change the parameter and return the hastings ratio.
     */
	public double doOperation() {
		
		double logq = 0;
		
		double currentValue;
		double newValue;
		
		final Bounds<Double> bounds = variable.getBounds();
		final int dim = variable.getSize();
		
		if (updateAllIndependently) {
			for (int i = 0; i < dim; i++) {
				
				//both current and new value of the variable needed for the hastings ratio
				currentValue = variable.getValue(i);
				
				newValue = gamma.nextGamma();
				while (newValue == 0.0) {
					newValue = gamma.nextGamma();
				}
				
				//System.err.println("newValue: " + newValue + " - logPdf: " + gamma.logPdf(newValue));
				
				logq += (gamma.logPdf(currentValue) - gamma.logPdf(newValue));
				
				if (newValue < bounds.getLowerLimit(i) || newValue > bounds.getUpperLimit(i)) {
                    throw new RuntimeException("proposed value outside boundaries");
                }
				
				variable.setValue(i, newValue);
				
			}
		}
		
		return logq;
	}
	
	public static dr.xml.XMLObjectParser PARSER = new dr.xml.AbstractXMLObjectParser() {
		
		public String getParserName() {
            return OPERATOR_NAME;
        }

		public Object parseXMLObject(XMLObject xo) throws XMLParseException {
			
			double weight = xo.getDoubleAttribute(WEIGHT);
			double shape = xo.getDoubleAttribute(SHAPE);
			double scale = xo.getDoubleAttribute(SCALE);
			
			if (! (shape > 0 && scale > 0)) {
				throw new XMLParseException("Shape and scale must be positive values.");
			}
			
			Parameter parameter = (Parameter) xo.getChild(Parameter.class);
			
			return new IndependentGammaSampler(parameter, new GammaDistribution(shape, scale), weight);
		}
		
		//************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

		public XMLSyntaxRule[] getSyntaxRules() {
			return rules;
		}
		
		private final XMLSyntaxRule[] rules = {
                AttributeRule.newDoubleRule(WEIGHT),
                AttributeRule.newDoubleRule(SHAPE),
                AttributeRule.newDoubleRule(SCALE),
                new ElementRule(Parameter.class)
        };

		public String getParserDescription() {
			return "This element returns an independence sampler from a provided gamma prior.";
		}

		public Class getReturnType() {
			return MCMCOperator.class;
		}
		
	};

}
