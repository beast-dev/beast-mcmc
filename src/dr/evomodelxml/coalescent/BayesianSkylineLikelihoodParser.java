/*
 * BayesianSkylineLikelihoodParser.java
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

package dr.evomodelxml.coalescent;

import dr.evomodel.coalescent.BayesianSkylineLikelihood;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.logging.Logger;

/**
 */
public class BayesianSkylineLikelihoodParser extends AbstractXMLObjectParser {

    public static final String SKYLINE_LIKELIHOOD = "generalizedSkyLineLikelihood";
    public static final String POPULATION_SIZES = "populationSizes";
    public static final String GROUP_SIZES = "groupSizes";

    public static final String TYPE = "type";
    public static final String STEPWISE = "stepwise";
    public static final String LINEAR = "linear";
    public static final String EXPONENTIAL = "exponential";

    public String getParserName() {
        return SKYLINE_LIKELIHOOD;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        XMLObject cxo = xo.getChild(POPULATION_SIZES);
        Parameter param = (Parameter) cxo.getChild(Parameter.class);

        cxo = xo.getChild(GROUP_SIZES);
        Parameter param2 = (Parameter) cxo.getChild(Parameter.class);

        cxo = xo.getChild(CoalescentLikelihoodParser.POPULATION_TREE);
        TreeModel treeModel = (TreeModel) cxo.getChild(TreeModel.class);

        int type = BayesianSkylineLikelihood.LINEAR_TYPE;
        String typeName = LINEAR;
        if (xo.hasAttribute(LINEAR) && !xo.getBooleanAttribute(LINEAR)) {
            type = BayesianSkylineLikelihood.STEPWISE_TYPE;
            typeName = STEPWISE;
        }

        if (xo.hasAttribute(TYPE)) {
            if (xo.getStringAttribute(TYPE).equalsIgnoreCase(STEPWISE)) {
                type = BayesianSkylineLikelihood.STEPWISE_TYPE;
                typeName = STEPWISE;
            } else if (xo.getStringAttribute(TYPE).equalsIgnoreCase(LINEAR)) {
                type = BayesianSkylineLikelihood.LINEAR_TYPE;
                typeName = LINEAR;
            } else if (xo.getStringAttribute(TYPE).equalsIgnoreCase(EXPONENTIAL)) {
                type = BayesianSkylineLikelihood.EXPONENTIAL_TYPE;
                typeName = EXPONENTIAL;
            } else throw new XMLParseException("Unknown Bayesian Skyline type: " + xo.getStringAttribute(TYPE));
        }

        if (param2.getDimension() > (treeModel.getExternalNodeCount()-1)) {
            throw new XMLParseException("There are more groups (" + param2.getDimension()
                    + ") than coalescent nodes in the tree (" + (treeModel.getExternalNodeCount()-1) + ").");
        }

        Logger.getLogger("dr.evomodel").info("Bayesian skyline plot: " + param.getDimension() + " " + typeName + " control points");

        return new BayesianSkylineLikelihood(treeModel, param, param2, type);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element represents the likelihood of the tree given the population size vector.";
    }

    public Class getReturnType() {
        return BayesianSkylineLikelihood.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules;{
        rules = new XMLSyntaxRule[]{
                new XORRule(
                        AttributeRule.newBooleanRule(LINEAR),
                        AttributeRule.newStringRule(TYPE)
                ),
                new ElementRule(POPULATION_SIZES, new XMLSyntaxRule[]{
                        new ElementRule(Parameter.class)
                }),
                new ElementRule(GROUP_SIZES, new XMLSyntaxRule[]{
                        new ElementRule(Parameter.class)
                }),
                new ElementRule(CoalescentLikelihoodParser.POPULATION_TREE, new XMLSyntaxRule[]{
                        new ElementRule(TreeModel.class)
                }),
        };
    }
}
