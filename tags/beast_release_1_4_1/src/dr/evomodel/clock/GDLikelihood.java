/*
 * GDLikelihood.java
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

package dr.evomodel.clock;

import dr.evomodel.tree.TreeModel;
import dr.math.GammaDistribution;
import dr.xml.*;
import dr.inference.model.Parameter;

/**
 * Calculates the likelihood of a set of rate changes in a tree, assuming a gamma distributed 
 * change in rate at each node, with a mean of the previous rate and a variance proportional to branch length.
 *
 * @author Alexei Drummond
 *
 * @version $Id: GDLikelihood.java,v 1.11 2005/05/24 20:25:57 rambaut Exp $
 */
public class GDLikelihood extends RateChangeLikelihood {
		
	public static final String GD_LIKELIHOOD = "GDLikelihood";	
	public static final String STDEV = "stdev";

	public GDLikelihood(TreeModel tree, Parameter ratesParameter, double stdev, int rootModel, boolean isEpisodic) {
		
		super("Gamma Distributed", tree, ratesParameter, rootModel, isEpisodic);
		this.unitVariance = stdev * stdev;
	}
	
    /**
     * @return the log likelihood of the rate change from the parent to the child.
     */
    double branchRateChangeLogLikelihood(double parentRate, double childRate, double time) {


        double variance = unitVariance;
        if (!isEpisodic()) {
            variance *= time;
        }

        double mean = parentRate;
        double shape = (mean * mean) / variance;
        double scale = variance / mean;

        return GammaDistribution.logPdf(childRate, shape, scale);
    }


	public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {
		
		public String getParserName() { return GD_LIKELIHOOD; }
		
		public Object parseXMLObject(XMLObject xo) throws XMLParseException {
			
			TreeModel tree = (TreeModel)xo.getChild(TreeModel.class);

            Parameter ratesParameter = (Parameter)xo.getSocketChild(RATES);

			double stdev = xo.getDoubleAttribute(STDEV);
            boolean episodic = xo.getBooleanAttribute(EPISODIC);

            String rootModelString = MEAN_OF_CHILDREN;
            int rootModel = ROOT_RATE_MEAN_OF_CHILDREN;
            if (xo.hasAttribute(ROOT_MODEL)) {
                rootModelString = xo.getStringAttribute(ROOT_MODEL);
                if (rootModelString.equals(MEAN_OF_CHILDREN)) rootModel = ROOT_RATE_MEAN_OF_CHILDREN;
                if (rootModelString.equals(MEAN_OF_ALL)) rootModel = ROOT_RATE_MEAN_OF_ALL;
                if (rootModelString.equals(EQUAL_TO_CHILD)) rootModel = ROOT_RATE_EQUAL_TO_CHILD;
                if (rootModelString.equals(IGNORE_ROOT)) rootModel = ROOT_RATE_IGNORE_ROOT;
                if (rootModelString.equals(NONE)) rootModel = ROOT_RATE_NONE;
            }

            System.out.println("Using auto-correlated relaxed clock model.");
            System.out.println("  parametric model = exponential distribution");
            System.out.println("  root rate model = " + rootModelString);

            return new GDLikelihood(tree, ratesParameter, stdev, rootModel, episodic);
		}
		
		//************************************************************************
		// AbstractXMLObjectParser implementation
		//************************************************************************
		
		public String getParserDescription() {
			return 
				"The likelihood of a set of rate changes in a tree, assuming "+
				"a gamma-distributed change in rate at each node, with a " + 
				"mean of the previous rate and a given variance (variance can be optionally proportional to " +
				"branch length).";
		}

		public Class getReturnType() { return GDLikelihood.class; }

		public XMLSyntaxRule[] getSyntaxRules() { return rules; }
	
		private XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
			new ElementRule(TreeModel.class),
            new ElementRule(RATES, Parameter.class, "The branch rates parameter", false),
			AttributeRule.newDoubleRule(STDEV, false, "The unit stdev of the model. The variance is scaled by the branch length to get the actual variance in the non-episodic version of the model."),
            AttributeRule.newStringRule(ROOT_MODEL, true, "specify the rate model to use at the root. Should be one of: 'meanOfChildren', 'meanOfAll', 'equalToChild', 'ignoreRoot' or 'none'."),
            AttributeRule.newBooleanRule(EPISODIC, false, "true if model is branch length independent, false if length-dependent.")
		};
	};
	
	double unitVariance = 1.0;
}