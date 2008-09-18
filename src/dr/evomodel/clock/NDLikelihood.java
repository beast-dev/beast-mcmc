/*
 * NDLikelihood.java
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
import dr.inference.model.Parameter;
import dr.math.distributions.NormalDistribution;
import dr.xml.*;


/**
 * Calculates the likelihood of a set of rate changes in a tree, assuming a normally distributed
 * change in rate at each node, with a mean of the previous rate and a variance proportional to branch length.
 *
 * @author Alexei Drummond
 * @version $Id: NDLikelihood.java,v 1.11 2005/05/24 20:25:57 rambaut Exp $
 */
public class NDLikelihood extends RateChangeLikelihood {

    public static final String ND_LIKELIHOOD = "NDLikelihood";
    public static final String VARIANCE = "variance";
    public static final String ROOTRATE = "rootRate";


    public NDLikelihood(TreeModel tree, Parameter ratesParameter, Parameter variance, Parameter rootRate, int rootModel, boolean isEpisodic, boolean isNormalized, boolean isLogSpace) {

        super((isLogSpace)?"LogNormally Distributed":"Normally Distributed", tree, ratesParameter, rootRate, rootModel, isEpisodic, isNormalized);

        this.isLogSpace = isLogSpace;
        this.variance = variance;

        addParameter(variance);
    }

    /**
     * @return the log likelihood of the rate change from the parent to the child.
     */
    double branchRateChangeLogLikelihood(double parentRate, double childRate, double time) {
        double var = variance.getParameterValue(0);

        if (isEpisodic())
            var *= time;

        if (isLogSpace){
            double logParentRate = Math.log(parentRate);
            double logChildRate = Math.log(childRate);

            return NormalDistribution.logPdf(logChildRate, logParentRate - (var /2.), Math.sqrt(var)) - logChildRate;

        } else{
            return NormalDistribution.logPdf(childRate, parentRate, Math.sqrt(var));
        }
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return ND_LIKELIHOOD;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            TreeModel tree = (TreeModel) xo.getChild(TreeModel.class);

            Parameter ratesParameter = (Parameter) xo.getElementFirstChild(RATES);

            Parameter variance = (Parameter)xo.getElementFirstChild(VARIANCE);

            Parameter rootRate = null;
             if (xo.hasChildNamed(ROOTRATE))
                 rootRate = (Parameter)xo.getElementFirstChild(ROOTRATE);

            boolean isEpisodic = xo.getBooleanAttribute(EPISODIC);

            boolean isNormalized = false;
            if (xo.hasAttribute(NORMALIZED))
                isNormalized = xo.getBooleanAttribute(NORMALIZED);

             boolean isLogSpace = false;
            if (xo.hasAttribute(LOGSPACE))
                isLogSpace = xo.getBooleanAttribute(LOGSPACE);

            String rootModelString = MEAN_OF_CHILDREN;
            int rootModel = ROOT_RATE_MEAN_OF_CHILDREN;
            if (xo.hasAttribute(ROOT_MODEL)) {
                rootModelString = xo.getStringAttribute(ROOT_MODEL);
                if (rootModelString.equals(MEAN_OF_CHILDREN)) rootModel = ROOT_RATE_MEAN_OF_CHILDREN;
                if (rootModelString.equals(MEAN_OF_ALL)) rootModel = ROOT_RATE_MEAN_OF_ALL;
                if (rootModelString.equals(EQUAL_TO_CHILD)) rootModel = ROOT_RATE_EQUAL_TO_CHILD;
                if (rootModelString.equals(IGNORE_ROOT)) rootModel = ROOT_RATE_IGNORE_ROOT;
                if (rootModelString.equals(PARAMETER)) rootModel = ROOT_RATE_PARAMETER;
                if (rootModelString.equals(NONE)) rootModel = ROOT_RATE_NONE;
            }

            return new NDLikelihood(tree, ratesParameter, variance,rootRate,  rootModel, isEpisodic, isNormalized, isLogSpace);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return
                    "This element returns an object that can calculate the likelihood " +
                            "of rate changes in a tree under the assumption of " +
                            "(log)normally distributed rate changes among lineages. " +
                            "Specifically, each branch is assumed to draw a rate from a " +
                            "(log)normal distribution with mean of the rate in the " +
                            "parent branch and the given standard deviation (the variance can be optionally proportional to " +
                            "branch length).";
        }

        public Class getReturnType() {
            return NDLikelihood.class;
        }


        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                new ElementRule(TreeModel.class),
                new ElementRule(RATES, Parameter.class, "The branch rates parameter", false),
                AttributeRule.newStringRule(ROOT_MODEL, true, "specify the rate model to use at the root. Should be one of: 'meanOfChildren', 'meanOfAll', 'equalToChild', 'ignoreRoot', 'parameterRoot' or 'none'."),
                AttributeRule.newBooleanRule(EPISODIC, false, "true if model is branch length independent, false if length-dependent."),
                AttributeRule.newBooleanRule(NORMALIZED, true, "true if relative rates"),
                AttributeRule.newBooleanRule(LOGSPACE, true, "true if model considers the log of the rates."),
                new ElementRule(ROOTRATE, Parameter.class, "The root rate parameter", true),
                new ElementRule(VARIANCE, Parameter.class, "The variance of the (log)normal distribution"),
        };
    };

   private Parameter variance;
    boolean isLogSpace = false;
}