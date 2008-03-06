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
import dr.math.NormalDistribution;
import dr.xml.*;
import dr.inference.model.Parameter;


/**
 * Calculates the likelihood of a set of rate changes in a tree, assuming a normally distributed
 * change in rate at each node, with a mean of the previous rate and a variance proportional to branch length.
 *
 * @author Alexei Drummond
 * @version $Id: NDLikelihood.java,v 1.11 2005/05/24 20:25:57 rambaut Exp $
 */
public class NDLikelihood extends RateChangeLikelihood {

    public static final String ND_LIKELIHOOD = "NDLikelihood";
    public static final String STDEV = "stdev";

    public NDLikelihood(TreeModel tree, Parameter ratesParameter, double stdev, int rootModel, boolean isEpisodic) {

        super("Normally Distributed", tree, ratesParameter, rootModel, isEpisodic, false);

        this.stdev = stdev;
    }

    /**
     * @return the log likelihood of the rate change from the parent to the child.
     */
    double branchRateChangeLogLikelihood(double parentRate, double childRate, double time) {
        if (isEpisodic()) {
            return NormalDistribution.logPdf(childRate, parentRate, stdev);
        } else {
            double scaledStdev = Math.sqrt(stdev * stdev * time);
            return NormalDistribution.logPdf(childRate, parentRate, scaledStdev);
        }
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return ND_LIKELIHOOD;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            TreeModel tree = (TreeModel) xo.getChild(TreeModel.class);

            Parameter ratesParameter = (Parameter) xo.getSocketChild(RATES);

            double stdev = xo.getDoubleAttribute(STDEV);
            boolean isEpisodic = xo.getBooleanAttribute(EPISODIC);

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

            NDLikelihood ed = new NDLikelihood(tree, ratesParameter, stdev, rootModel, isEpisodic);

            return ed;
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return
                    "This element returns an object that can calculate the likelihood " +
                            "of rate changes in a tree under the assumption of " +
                            "normally distributed rate changes among lineages. " +
                            "Specifically, each branch is assumed to draw a rate from a " +
                            "normal distribution with mean of the rate in the " +
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
                AttributeRule.newDoubleRule(STDEV, false, "The unit stdev of the model. The variance is scaled by the branch length to get the actual variance in the non-episodic version of the model."),
                AttributeRule.newStringRule(ROOT_MODEL, true, "specify the rate model to use at the root. Should be one of: 'meanOfChildren', 'meanOfAll', 'equalToChild', 'ignoreRoot' or 'none'."),
                AttributeRule.newBooleanRule(EPISODIC, false, "true if model is branch length independent, false if length-dependent.")
        };
    };

    double stdev = 1.0;
}