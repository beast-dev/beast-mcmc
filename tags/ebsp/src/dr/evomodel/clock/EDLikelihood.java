/*
 * EDLikelihood.java
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
import dr.math.ExponentialDistribution;
import dr.xml.*;
import dr.inference.model.Parameter;


/**
 * Calculates the likelihood of a set of rate changes in a tree, assuming an exponentially distributed
 * change in rate at each node, with a mean of the previous rate. This model ignores the branch lengths.
 *
 * @author Alexei Drummond
 * @version $Id: EDLikelihood.java,v 1.14 2005/05/24 20:25:57 rambaut Exp $
 */
public class EDLikelihood extends RateChangeLikelihood {

    public static final String ED_LIKELIHOOD = "EDLikelihood";

    public EDLikelihood(TreeModel tree, Parameter ratesParameter, int rootRatePrior, boolean isNormalized) {

        super("Exponentially Distributed", tree, ratesParameter, rootRatePrior, true, isNormalized);


    }

    /**
     * @return the log likelihood of the rate change from the parent to the child.
     */
    double branchRateChangeLogLikelihood(double parentRate, double childRate, double time) {
        return ExponentialDistribution.logPdf(childRate, 1.0 / parentRate);
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return ED_LIKELIHOOD;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            TreeModel tree = (TreeModel) xo.getChild(TreeModel.class);

            Parameter ratesParameter = (Parameter) xo.getSocketChild(RATES);

            boolean isNormalized = xo.getBooleanAttribute(NORMALIZED);

            String rootModelString = MEAN_OF_CHILDREN;
            int rootModel = ROOT_RATE_MEAN_OF_CHILDREN;
            if (xo.hasAttribute(ROOT_MODEL)) {
                rootModelString = xo.getStringAttribute(ROOT_MODEL);
                if (rootModelString.equals(MEAN_OF_CHILDREN)) rootModel = ROOT_RATE_MEAN_OF_CHILDREN;
                if (rootModelString.equals(MEAN_OF_ALL)) rootModel = ROOT_RATE_MEAN_OF_ALL;
                if (rootModelString.equals(EQUAL_TO_CHILD)) rootModel = ROOT_RATE_EQUAL_TO_CHILD;
                if (rootModelString.equals(IGNORE_ROOT)) rootModel = ROOT_RATE_IGNORE_ROOT;
                if (rootModelString.equals(FIXED_ROOT)) rootModel = ROOT_RATE_FIXED_ROOT;
                if (rootModelString.equals(NONE)) rootModel = ROOT_RATE_NONE;
            }

            System.out.println("Using auto-correlated relaxed clock model.");
            System.out.println("  parametric model = exponential distribution");
            System.out.println("  root rate model = " + rootModelString);

            return new EDLikelihood(tree, ratesParameter, rootModel, isNormalized);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return
                    "This element returns an object that can calculate the likelihood " +
                            "of rate changes in a tree under the assumption of " +
                            "exponentially distributed rate changes among lineages. " +
                            "Specifically, each branch is assumed to draw a rate from an " +
                            "exponential distribution with mean of the rate in the " +
                            "parent branch.";
        }

        public Class getReturnType() {
            return EDLikelihood.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                new ElementRule(TreeModel.class),
                AttributeRule.newStringRule(ROOT_MODEL, true, "specify the rate model to use at the root. Should be one of: 'meanOfChildren', 'meanOfAll', 'equalToChild', 'ignoreRoot', 'fixedRoot' or 'none'."),
                AttributeRule.newBooleanRule(NORMALIZED, false, "true if relative rates"),
                new ElementRule(RATES, Parameter.class, "The branch rates parameter", false)
        };
    };


}