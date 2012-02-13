/*
 * GaussianProcessSkytrackLikelihoodParser.java
 *
 * Copyright (c) 2002-2011 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

import dr.evolution.tree.Tree;
//import dr.evomodel.coalescent.GMRFSkyrideLikelihood;
import dr.evomodel.coalescent.GaussianProcessSkytrackLikelihood;
import dr.evomodel.tree.TreeModel;
//import dr.inference.model.MatrixParameter;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;

/**
 *
 */
public class GaussianProcessSkytrackLikelihoodParser extends AbstractXMLObjectParser {

//	public static final String SKYLINE_LIKELIHOOD = "gmrfSkyrideLikelihood";
    public static final String SKYTRACK_LIKELIHOOD = "gpSkytrackLikelihood";

    public static final String LAMBDA_BOUND_PARAMETER = "lambdaBoundParameter";
	public static final String POPULATION_PARAMETER = "populationSizes";
    public static final String CHANGE_POINTS="changePoints";
//	public static final String GROUP_SIZES = "groupSizes";
	public static final String PRECISION_PARAMETER = "precisionParameter";
	public static final String POPULATION_TREE = "populationTree";
	public static final String LAMBDA_PARAMETER = "lambdaParameter";
	public static final String BETA_PARAMETER = "betaParameter";
    public static final String ALPHA_PARAMETER = "alphaParameter";
//	public static final String COVARIATE_MATRIX = "covariateMatrix";
	public static final String RANDOMIZE_TREE = "randomizeTree";
//	public static final String TIME_AWARE_SMOOTHING = "timeAwareSmoothing";

    public static final String RESCALE_BY_ROOT_ISSUE = "rescaleByRootHeight";
//    public static final String GRID_POINTS = "gridPoints";
//    public static final String OLD_SKYRIDE = "oldSkyride";      //True=No multiple loci
    public static final String NUM_GRID_POINTS = "numGridPoints";
//    public static final String CUT_OFF = "cutOff";
//    public static final String PHI_PARAMETER = "phiParameter";

//    public static final String LATENT_PARAMETER = "latentPointParameter";
    

    public String getParserName() {
       return SKYTRACK_LIKELIHOOD;
    }



    public Object parseXMLObject(XMLObject xo) throws XMLParseException {


        XMLObject cxo = xo.getChild(POPULATION_PARAMETER);
        Parameter popParameter = (Parameter) cxo.getChild(Parameter.class);

        cxo = xo.getChild(PRECISION_PARAMETER);
        Parameter precParameter = (Parameter) cxo.getChild(Parameter.class);

//        cxo = xo.getChild(LAMBDA_BOUND_PARAMETER);
//        Parameter lambda_bound = (Parameter) cxo.getChild(Parameter.class);
//
//        cxo = xo.getChild(LAMBDA_PARAMETER);
//        Parameter lambda_parameter = (Parameter) cxo.getChild(Parameter.class);


        cxo = xo.getChild(POPULATION_TREE);

        List<Tree> treeList = new ArrayList<Tree>();
        for (int i = 0; i < cxo.getChildCount(); i++) {
            Object testObject = cxo.getChild(i);
            if (testObject instanceof Tree) {
                treeList.add((TreeModel) testObject);
            }
        }

//        TreeModel treeModel = (TreeModel) cxo.getChild(TreeModel.class);

//        cxo = xo.getChild(GROUP_SIZES);
//        Parameter groupParameter = null;
//        if (cxo != null) {
//            groupParameter = (Parameter) cxo.getChild(Parameter.class);
//
//            if (popParameter.getDimension() != groupParameter.getDimension())
//                throw new XMLParseException("Population and group size parameters must have the same length");
//        }

        Parameter lambda_parameter;
        if (xo.getChild(LAMBDA_PARAMETER) != null) {
            cxo = xo.getChild(LAMBDA_PARAMETER);
            lambda_parameter = (Parameter) cxo.getChild(Parameter.class);
        } else {
            lambda_parameter = new Parameter.Default(1.0);
        }

        Parameter lambda_bound;
        if (xo.getChild(LAMBDA_BOUND_PARAMETER) != null) {
            cxo = xo.getChild(LAMBDA_BOUND_PARAMETER);
            lambda_bound = (Parameter) cxo.getChild(Parameter.class);
        } else {
            lambda_bound = new Parameter.Default(1.0);
        }

        Parameter alpha_parameter;
        if (xo.getChild(ALPHA_PARAMETER) != null) {
            cxo = xo.getChild(ALPHA_PARAMETER);
            alpha_parameter = (Parameter) cxo.getChild(Parameter.class);
        } else {
            alpha_parameter = new Parameter.Default(0.001);
        }

        Parameter beta_parameter;
               if (xo.getChild(BETA_PARAMETER) != null) {
                   cxo = xo.getChild(BETA_PARAMETER);
                   beta_parameter = (Parameter) cxo.getChild(Parameter.class);
               } else {
                   beta_parameter = new Parameter.Default(0.001);
               }

        Parameter change_points;
               if (xo.getChild(CHANGE_POINTS) != null) {
                   cxo = xo.getChild(CHANGE_POINTS);
                   change_points = (Parameter) cxo.getChild(Parameter.class);
               } else {
                   change_points = new Parameter.Default(0,1);
               }

                /*


        Parameter gridPoints = null;
        if (xo.getChild(GRID_POINTS) != null) {
            cxo = xo.getChild(GRID_POINTS);
            gridPoints = (Parameter) cxo.getChild(Parameter.class);
        }
        */


        Parameter numGridPoints = null;
        if (xo.getChild(NUM_GRID_POINTS) != null) {
            cxo = xo.getChild(NUM_GRID_POINTS);
            numGridPoints = (Parameter) cxo.getChild(Parameter.class);
        }

//        Parameter cutOff = null;
//        if (xo.getChild(CUT_OFF) != null) {
//            cxo = xo.getChild(CUT_OFF);
//            cutOff = (Parameter) cxo.getChild(Parameter.class);
//        }
//
//        Parameter phi = null;
//        if (xo.getChild(PHI_PARAMETER) != null) {
//            cxo = xo.getChild(PHI_PARAMETER);
//            phi = (Parameter) cxo.getChild(Parameter.class);
//        }
//
//
//        Parameter beta = null;
//        if (xo.getChild(BETA_PARAMETER) != null) {
//            cxo = xo.getChild(BETA_PARAMETER);
//            beta = (Parameter) cxo.getChild(Parameter.class);
//        }

//        MatrixParameter dMatrix = null;
//        if (xo.getChild(COVARIATE_MATRIX) != null) {
//            cxo = xo.getChild(COVARIATE_MATRIX);
//            dMatrix = (MatrixParameter) cxo.getChild(MatrixParameter.class);
//        }

//        boolean timeAwareSmoothing = GMRFSkyrideLikelihood.TIME_AWARE_IS_ON_BY_DEFAULT;
//        if (xo.hasAttribute(TIME_AWARE_SMOOTHING)) {
//            timeAwareSmoothing = xo.getBooleanAttribute(TIME_AWARE_SMOOTHING);
//        }

//        if ((dMatrix != null && beta == null) || (dMatrix == null && beta != null))
//            throw new XMLParseException("Must specify both a set of regression coefficients and a design matrix.");
//
//        if (dMatrix != null) {
//            if (dMatrix.getRowDimension() != popParameter.getDimension())
//                throw new XMLParseException("Design matrix row dimension must equal the population parameter length.");
//            if (dMatrix.getColumnDimension() != beta.getDimension())
//                throw new XMLParseException("Design matrix column dimension must equal the regression coefficient length.");
//        }

        if (xo.getAttribute(RANDOMIZE_TREE, false)) {
            for (Tree tree : treeList) {
                if (tree instanceof TreeModel) {
                    GaussianProcessSkytrackLikelihood.checkTree((TreeModel) tree);
                } else {
                    throw new XMLParseException("Can not randomize a fixed tree");
                }
            }
        }


//        XMLObject latentChild = xo.getChild(LATENT_PARAMETER);
//        Parameter latentPoints = (Parameter) latentChild.getChild(Parameter.class);


        boolean rescaleByRootHeight = xo.getAttribute(RESCALE_BY_ROOT_ISSUE, true);

//        Logger.getLogger("dr.evomodel").info("The " + SKYTRACK_LIKELIHOOD + " has " +
//                (timeAwareSmoothing ? "time aware smoothing" : "uniform smoothing"));


             return new GaussianProcessSkytrackLikelihood(treeList, precParameter,
                 rescaleByRootHeight, numGridPoints, lambda_bound, lambda_parameter, popParameter,alpha_parameter,beta_parameter, change_points);

    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element represents the likelihood of the tree given the population size vector.";
    }

    public Class getReturnType() {
        return GaussianProcessSkytrackLikelihood.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(POPULATION_PARAMETER, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }),
            new ElementRule(CHANGE_POINTS, new XMLSyntaxRule[]{
                              new ElementRule(Parameter.class)
            }),
            new ElementRule(PRECISION_PARAMETER, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }),
            new ElementRule(NUM_GRID_POINTS, new XMLSyntaxRule[]{
                               new ElementRule(Parameter.class)
                       }),
            new ElementRule(LAMBDA_BOUND_PARAMETER, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }),
            new ElementRule(LAMBDA_PARAMETER, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }),
            new ElementRule(ALPHA_PARAMETER, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }),
            new ElementRule(BETA_PARAMETER, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }),
//            new ElementRule(PHI_PARAMETER, new XMLSyntaxRule[]{
//                    new ElementRule(Parameter.class)
//            }, true), // Optional
            new ElementRule(POPULATION_TREE, new XMLSyntaxRule[]{
                    new ElementRule(TreeModel.class, 1, Integer.MAX_VALUE)
            }),
//            new ElementRule(GROUP_SIZES, new XMLSyntaxRule[]{
//                    new ElementRule(Parameter.class)
//            }, true),
//            new ElementRule(LATENT_PARAMETER,Parameter.class),

            AttributeRule.newBooleanRule(RESCALE_BY_ROOT_ISSUE, true),
            AttributeRule.newBooleanRule(RANDOMIZE_TREE, true),
//            AttributeRule.newBooleanRule(TIME_AWARE_SMOOTHING, true),
//            AttributeRule.newBooleanRule(OLD_SKYRIDE, true)
    };

}
