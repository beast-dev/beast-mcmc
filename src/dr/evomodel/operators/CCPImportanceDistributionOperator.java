/*
 * CCPImportanceDistributionOperator.java
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

/**
 *
 */
package dr.evomodel.operators;

import dr.evomodel.tree.ConditionalCladeFrequency;
import dr.evomodel.tree.TreeModel;
import dr.xml.*;

/**
 * @author Sebastian Hoehna
 */
// Cleaning out untouched stuff. Can be resurrected if needed
@Deprecated
public class CCPImportanceDistributionOperator extends
        AbstractImportanceDistributionOperator {

    public static final String CCP_IMPORTANCE_DISTRIBUTION_OPERATOR = "CCPImportanceDistributionOperator";

    /**
     *
     */
    public CCPImportanceDistributionOperator(TreeModel tree, double weight,
                                             int samples, int sampleEvery, double epsilon) {
        super(tree, weight, samples, sampleEvery);

        probabilityEstimater = new ConditionalCladeFrequency(tree, epsilon);
    }

    /**
     *
     */
    public CCPImportanceDistributionOperator(TreeModel tree, double weight) {
        super(tree, weight);
        double epsilon = 1 - Math.pow(0.5, 1.0 / 10000);
        probabilityEstimater = new ConditionalCladeFrequency(tree, epsilon);
    }

    /*
      * (non-Javadoc)
      *
      * @see dr.inference.operators.AbstractImportanceSampler#getOperatorName()
      */
    @Override
    public String getOperatorName() {
        return CCP_IMPORTANCE_DISTRIBUTION_OPERATOR;
    }

    /*
      * (non-Javadoc)
      *
      * @see
      * dr.inference.operators.AbstractImportanceSampler#getPerformanceSuggestion
      * ()
      */
    @Override
    public String getPerformanceSuggestion() {
        if (getAcceptanceProbability() < getMinimumGoodAcceptanceLevel()) {
            return "Try to increase the sample size and/or the steps between each sample.";
        }
        return "";
    }


    // Sebastian
//	public void printClades(){
//		probabilityEstimater.printClades();
//	}
    // Sebastian

    public static XMLObjectParser CCP_IMPORTANCE_DISTRIBUTION_OPERATOR_PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return CCP_IMPORTANCE_DISTRIBUTION_OPERATOR;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
            double weight = xo.getDoubleAttribute("weight");
            int samples = xo.getIntegerAttribute("samples");

            double epsilon = 1 - Math.pow(0.5, 1.0 / samples);
            if (xo.hasAttribute("epsilon")) {
                epsilon = xo.getDoubleAttribute("epsilon");
            }

            int sampleEvery = 10;
            if (xo.hasAttribute("sampleEvery")) {
                sampleEvery = xo.getIntegerAttribute("sampleEvery");
            }

            return new CCPImportanceDistributionOperator(treeModel, weight,
                    samples, sampleEvery, epsilon);
        }

        //**********************************************************************
        // **
        // AbstractXMLObjectParser implementation
        //**********************************************************************
        // **

        public String getParserDescription() {
            return "This element represents an operator proposing trees from an importance distribution which is created by the conditional clade probabilities.";
        }

        public Class getReturnType() {
            return CCPImportanceDistributionOperator.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                AttributeRule.newDoubleRule("weight"),
                AttributeRule.newIntegerRule("samples"),
                AttributeRule.newIntegerRule("sampleEvery", true),
                AttributeRule.newDoubleRule("epsilon", true),
                new ElementRule(TreeModel.class)};

    };
}
