/*
 * GMRFSkyrideLikelihoodParser.java
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

import dr.evolution.coalescent.IntervalList;
import dr.evomodel.coalescent.*;
import dr.evolution.tree.Tree;
import dr.evomodel.tree.TreeModel;
import dr.inference.glm.GeneralizedLinearModel;
import dr.inference.model.MatrixParameter;
import dr.inference.model.DesignMatrix;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;

/**
 *
 */
public class GMRFSkyrideLikelihoodParser extends AbstractXMLObjectParser {

    public static final String SKYLINE_LIKELIHOOD = "gmrfSkyrideLikelihood";
    public static final String SKYRIDE_LIKELIHOOD = "skyrideLikelihood";
    public static final String SKYGRID_LIKELIHOOD = "gmrfSkyGridLikelihood";

    public static final String POPULATION_PARAMETER = "populationSizes";
    public static final String GROUP_SIZES = "groupSizes";
    public static final String PRECISION_PARAMETER = "precisionParameter";
    public static final String POPULATION_TREE = "populationTree";
    public static final String INTERVALS = "intervals";
    public static final String BUILD_MAPPING = "intervalNodeMapping";
    public static final String LAMBDA_PARAMETER = "lambdaParameter";
    public static final String BETA_PARAMETER = "betaParameter";
    public static final String DELTA_PARAMETER = "deltaParameter";
    public static final String SINGLE_BETA = "singleBeta";
    public static final String COVARIATE_MATRIX = "covariateMatrix";
    public static final String RANDOMIZE_TREE = "randomizeTree";
    public static final String TIME_AWARE_SMOOTHING = "timeAwareSmoothing";

    public static final String RESCALE_BY_ROOT_ISSUE = "rescaleByRootHeight";
    public static final String GRID_POINTS = "gridPoints";
    public static final String OLD_SKYRIDE = "oldSkyride";
    public static final String NUM_GRID_POINTS = "numGridPoints";
    public static final String CUT_OFF = "cutOff";
    public static final String PHI_PARAMETER = "phiParameter";
    public static final String PLOIDY = "ploidy";
    public static final String COVARIATES = "covariates";
    public static final String COLUMN_MAJOR = "columnMajor";
    public static final String FIRST_OBSERVED_INDEX = "firstObservedIndex";
    public static final String LAST_OBSERVED_INDEX = "lastObservedIndex";
    public static final String COV_PREC_PARAM = "covariatePrecision";
    public static final String COV_PREC_REC = "covariatePrecisionRecent";
    public static final String COV_PREC_DIST = "covariatePrecisionDistant";
    public static final String REC_INDICES = "covIndicesMissingRecent";
    public static final String DIST_INDICES = "covIndicesMissingDistant";
    public static final String GLM_MODEL = "glmModel";
    public static final String USE_GLM_MODEL = "useGlmModel";

    public String getParserName() {
        return SKYLINE_LIKELIHOOD;
    }

    public String[] getParserNames() {
        return new String[]{getParserName(), SKYRIDE_LIKELIHOOD, SKYGRID_LIKELIHOOD}; // cannot duplicate
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        XMLObject cxo = xo.getChild(POPULATION_PARAMETER);
        Parameter popParameter = (Parameter) cxo.getChild(Parameter.class);

        cxo = xo.getChild(PRECISION_PARAMETER);
        Parameter precParameter = (Parameter) cxo.getChild(Parameter.class);

        boolean buildIntervalNodeMapping = xo.getAttribute(BUILD_MAPPING, false);

        List<IntervalList> intervalsList = new ArrayList<IntervalList>();

        List<Tree> treeList = new ArrayList<Tree>();
        if(xo.getChild(POPULATION_TREE) != null) {
            cxo = xo.getChild(POPULATION_TREE);
            for (int i = 0; i < cxo.getChildCount(); i++){
                Object testObject = cxo.getChild(i);
                if (testObject instanceof Tree) {
                    treeList.add((TreeModel) testObject);

//                    TreeIntervals treeIntervals;
//                    try {
//                        treeIntervals = new TreeIntervals((Tree) testObject, null, null);
//                    } catch (TreeUtils.MissingTaxonException mte) {
//                        throw new XMLParseException("Taxon, " + mte + ", in " + getParserName() + " was not found in the tree.");
//                    }
//                    intervalsList.add(treeIntervals);
                }
            }
        }

        /*
        List<Tree> treeList = new ArrayList<Tree>();
        for (int i = 0; i < cxo.getChildCount(); i++) {
            Object testObject = cxo.getChild(i);
            if (testObject instanceof Tree) {
                treeList.add((TreeModel) testObject);
            }
        } */

        if (xo.getChild(INTERVALS) != null) {
            cxo = xo.getChild(INTERVALS);
            intervalsList = new ArrayList<IntervalList>();
            for (int i = 0; i < cxo.getChildCount(); i++) {
                Object testObject = cxo.getChild(i);
                if (testObject instanceof IntervalList) {
                    intervalsList.add((IntervalList) testObject);
                }
            }
        }

        cxo = xo.getChild(GROUP_SIZES);
        Parameter groupParameter = null;
        if (cxo != null) {
            groupParameter = (Parameter) cxo.getChild(Parameter.class);

            if (popParameter.getDimension() != groupParameter.getDimension())
                throw new XMLParseException("Population and group size parameters must have the same length");
        }

        Parameter lambda;
        if (xo.getChild(LAMBDA_PARAMETER) != null) {
            cxo = xo.getChild(LAMBDA_PARAMETER);
            lambda = (Parameter) cxo.getChild(Parameter.class);
        } else {
            lambda = new Parameter.Default(LAMBDA_PARAMETER, 1.0);
        }

        Parameter gridPoints = null;
        if (xo.getChild(GRID_POINTS) != null) {
            cxo = xo.getChild(GRID_POINTS);
            gridPoints = (Parameter) cxo.getChild(Parameter.class);
        }

        Parameter numGridPoints = null;
        if (xo.getChild(NUM_GRID_POINTS) != null) {
            cxo = xo.getChild(NUM_GRID_POINTS);
            numGridPoints = (Parameter) cxo.getChild(Parameter.class);
        }

        Parameter cutOff = null;
        if (xo.getChild(CUT_OFF) != null) {
            cxo = xo.getChild(CUT_OFF);
            cutOff = (Parameter) cxo.getChild(Parameter.class);
        }

        Parameter phi = null;
        if (xo.getChild(PHI_PARAMETER) != null) {
            cxo = xo.getChild(PHI_PARAMETER);
            phi = (Parameter) cxo.getChild(Parameter.class);
        }

        List<Parameter> firstObservedIndex = null;
        if (xo.hasChildNamed(FIRST_OBSERVED_INDEX)) {
            firstObservedIndex = new ArrayList<Parameter>();
            cxo = xo.getChild(FIRST_OBSERVED_INDEX);
            final int numInd = cxo.getChildCount();

            for(int i=0; i< numInd; ++i) {
                firstObservedIndex.add((Parameter) cxo.getChild(i));
            }
        }

        List<Parameter> lastObservedIndex = null;
        if (xo.hasChildNamed(LAST_OBSERVED_INDEX)) {
            lastObservedIndex = new ArrayList<Parameter>();
            cxo = xo.getChild(LAST_OBSERVED_INDEX);
            final int numObsInd = cxo.getChildCount();

            for(int i=0; i< numObsInd; ++i) {
                lastObservedIndex.add((Parameter) cxo.getChild(i));
            }
        }

        Parameter ploidyFactors = null;
        if (xo.getChild(PLOIDY) != null) {
            cxo = xo.getChild(PLOIDY);
            ploidyFactors = (Parameter) cxo.getChild(Parameter.class);
        } else {
            if (intervalsList.size() != 0) {
                ploidyFactors = new Parameter.Default(PLOIDY, intervalsList.size());
                for(int i = 0; i < intervalsList.size(); i++){
                    ploidyFactors.setParameterValue(i, 1.0);
                }
            } else {
                ploidyFactors = new Parameter.Default(PLOIDY, treeList.size());
                for (int i = 0; i < treeList.size(); i++) {
                    ploidyFactors.setParameterValue(i, 1.0);
                }
            }
        }

        Parameter betaParameter = null;
        if (xo.hasChildNamed(SINGLE_BETA)) {
            betaParameter = (Parameter) xo.getElementFirstChild(SINGLE_BETA);
        }

        List<Parameter> betaList = null;
        if (xo.getChild(BETA_PARAMETER) != null) {
            betaList = new ArrayList<Parameter>();
            cxo = xo.getChild(BETA_PARAMETER);
            final int numBeta = cxo.getChildCount();
            for (int i = 0; i < numBeta; ++i) {
                betaList.add((Parameter) cxo.getChild(i));
            }
        }

        List<Parameter> deltaList = new ArrayList<Parameter>();
        if (xo.getChild(DELTA_PARAMETER) != null) {
            cxo = xo.getChild(DELTA_PARAMETER);
            final int numDelta = cxo.getChildCount();
            if(numDelta != betaList.size()){
                throw new XMLParseException("Cannot have different number of delta and beta parameters");
            }
            for (int i = 0; i < numDelta; ++i) {
                deltaList.add((Parameter) cxo.getChild(i));
            }
        }else{
            deltaList = null;
        }

        MatrixParameter dMatrix = null;
        if (xo.getChild(COVARIATE_MATRIX) != null) {
            cxo = xo.getChild(COVARIATE_MATRIX);
            dMatrix = (MatrixParameter) cxo.getChild(MatrixParameter.class);
        }

        boolean timeAwareSmoothing = GMRFSkyrideLikelihood.TIME_AWARE_IS_ON_BY_DEFAULT;
        if (xo.hasAttribute(TIME_AWARE_SMOOTHING)) {
            timeAwareSmoothing = xo.getBooleanAttribute(TIME_AWARE_SMOOTHING);
        }

        if (dMatrix != null) {
            if (dMatrix.getRowDimension() != popParameter.getDimension())
                throw new XMLParseException("Design matrix row dimension must equal the population parameter length.");
            if (dMatrix.getColumnDimension() != betaParameter.getDimension())
                throw new XMLParseException("Design matrix column dimension must equal the regression coefficient length.");
        }

        List<Parameter> covPrecParamRecent = null;
        List<Parameter> covPrecParamDistant = null;

        if(xo.hasChildNamed(COV_PREC_REC)){
            covPrecParamRecent = new ArrayList<Parameter>();
            cxo = xo.getChild(COV_PREC_REC);
            for(int i = 0; i < cxo.getChildCount(); ++i){
                covPrecParamRecent.add((Parameter) cxo.getChild(i));
            }
        }

        if(xo.hasChildNamed(COV_PREC_DIST)){
            covPrecParamDistant = new ArrayList<Parameter>();
            cxo = xo.getChild(COV_PREC_DIST);
            for(int i = 0; i < cxo.getChildCount(); ++i){
                covPrecParamDistant.add((Parameter) cxo.getChild(i));
            }
        }

        if (xo.hasChildNamed(COV_PREC_PARAM)){
            if(firstObservedIndex != null) {
                covPrecParamRecent = new ArrayList<Parameter>();
            }
            if(lastObservedIndex != null) {
                covPrecParamDistant = new ArrayList<Parameter>();
            }
            cxo = xo.getChild(COV_PREC_PARAM);

            for(int i=0; i < cxo.getChildCount(); ++i){
                if(firstObservedIndex != null) {
                    covPrecParamRecent.add((Parameter) cxo.getChild(i));
                }
                if(lastObservedIndex != null) {
                    covPrecParamDistant.add((Parameter) cxo.getChild(i));
                }
            }
        }

        if((covPrecParamDistant == null && lastObservedIndex != null) || (covPrecParamDistant != null && lastObservedIndex == null)){
            throw new XMLParseException("Must specify both lastObservedIndex and covariatePrecision");
        }

        if((covPrecParamRecent == null && firstObservedIndex != null) || (covPrecParamRecent != null && firstObservedIndex == null)){
            throw new XMLParseException("Must specify both firstObservedIndex and covariatePrecision");
        }

        Parameter recentIndices = null;
        if (xo.getChild(REC_INDICES) != null) {
            cxo = xo.getChild(REC_INDICES);
            recentIndices = (Parameter) cxo.getChild(Parameter.class);
        }

        if(firstObservedIndex == null && recentIndices != null){
            throw new XMLParseException("Cannot specify covIndicesMissingRecent without specifying firstObservedIndex");
        }

        Parameter distantIndices = null;
        if (xo.getChild(DIST_INDICES) != null) {
            cxo = xo.getChild(DIST_INDICES);
            distantIndices = (Parameter) cxo.getChild(Parameter.class);
        }

        if(lastObservedIndex == null && distantIndices != null){
            throw new XMLParseException("Cannot specify covIndicesMissingDistant without specifying lastObservedIndex");
        }

        List<MatrixParameter> covariates = null;
        if (xo.hasChildNamed(COVARIATES)){
            covariates = new ArrayList<MatrixParameter>();
            cxo = xo.getChild(COVARIATES);
            final int numCov = cxo.getChildCount();

            for (int i = 0; i < numCov; ++i) {
                covariates.add((MatrixParameter) cxo.getChild(i));
            }
        }

        if ((covariates != null && betaList == null) ||
                (covariates == null &&  betaList != null))
            throw new XMLParseException("Must specify both a set of regression coefficients and a design matrix.");

        boolean useGlmModel = xo.getAttribute(USE_GLM_MODEL, false);

        if(useGlmModel) {

            GeneralizedLinearModel glm = (GeneralizedLinearModel) xo.getChild(GeneralizedLinearModel.class);
            covariates = new ArrayList<MatrixParameter>();
            betaList = new ArrayList<Parameter>();
            List<DesignMatrix> designMat = glm.getDesignMatrix();
            List<Parameter> indepParam = glm.getIndependentParameter();
            List<Parameter> indepParamDelta = glm.getIndependentParameterDelta();
            deltaList = new ArrayList<Parameter>();

            for(int i = 0; i < indepParam.get(0).getSize(); i++){
                MatrixParameter matParam = new MatrixParameter("covariate values", 1,designMat.get(0).getRowDimension());

                for(int j = 0; j < matParam.getRowDimension(); j++){
                    matParam.setParameterValue(0, j, designMat.get(0).getParameterValue(0, j));
                }
                covariates.add(matParam);

                Parameter betaParam = new Parameter.Default(1);
                betaParam.setParameterValue(0, indepParam.get(0).getParameterValue(i));
                betaList.add(betaParam);

                if(indepParamDelta != null){
                    Parameter deltaParam = new Parameter.Default(1);
                    deltaParam.setParameterValue(0, indepParamDelta.get(0).getParameterValue(i));
                    deltaList.add(deltaParam);
                }
            }

        }

        /*
        if (xo.getAttribute(RANDOMIZE_TREE, false)) {
            for (Tree tree : treeList) {
                if (tree instanceof TreeModel) {
                    GMRFSkyrideLikelihood.checkTree((TreeModel) tree);
                } else {
                    throw new XMLParseException("Can not randomize a fixed tree");
                }
            }
        }*/

        boolean rescaleByRootHeight = xo.getAttribute(RESCALE_BY_ROOT_ISSUE, true);

        Logger.getLogger("dr.evomodel").info("The " + SKYLINE_LIKELIHOOD + " has " +
                (timeAwareSmoothing ? "time aware smoothing" : "uniform smoothing"));

        if (xo.getAttribute(OLD_SKYRIDE, true) && xo.getName().compareTo(SKYGRID_LIKELIHOOD) != 0) {
            return new OldGMRFSkyrideLikelihood(treeList, popParameter, groupParameter, precParameter,
                    lambda, betaParameter, dMatrix, timeAwareSmoothing, rescaleByRootHeight, buildIntervalNodeMapping);

        } else {
            if (intervalsList.size() > 0) {
                if (xo.getChild(GRID_POINTS) != null) {
                    return new GMRFSkygridLikelihood(intervalsList, popParameter, groupParameter, precParameter,
                            lambda, betaParameter, dMatrix, timeAwareSmoothing, gridPoints, covariates, ploidyFactors,
                            firstObservedIndex, lastObservedIndex, covPrecParamRecent, covPrecParamDistant, recentIndices, distantIndices, betaList);
                } else {
                    return new GMRFSkygridLikelihood(intervalsList, popParameter, groupParameter, precParameter,
                            lambda, betaParameter, dMatrix, timeAwareSmoothing, cutOff.getParameterValue(0), (int) numGridPoints.getParameterValue(0), phi, ploidyFactors);
                }

            } else {
                if (xo.getChild(GRID_POINTS) != null) {
                    return new GMRFMultilocusSkyrideLikelihood(treeList, popParameter, groupParameter, precParameter,
                            lambda, betaParameter, dMatrix, timeAwareSmoothing, gridPoints, covariates, ploidyFactors,
                            firstObservedIndex, lastObservedIndex, covPrecParamRecent, covPrecParamDistant, recentIndices, distantIndices, betaList, deltaList);
                } else {
                    return new GMRFMultilocusSkyrideLikelihood(treeList, popParameter, groupParameter, precParameter,
                            lambda, betaParameter, dMatrix, timeAwareSmoothing, cutOff.getParameterValue(0), (int) numGridPoints.getParameterValue(0), phi, ploidyFactors);
                }
            }
        }
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element represents the likelihood of the tree given the population size vector.";
    }

    public Class getReturnType() {
        return GMRFSkyrideLikelihood.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(POPULATION_PARAMETER, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }),
            new ElementRule(PRECISION_PARAMETER, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }),
            new ElementRule(PHI_PARAMETER, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }, true), // Optional
            new OrRule(
                new ElementRule(INTERVALS, new XMLSyntaxRule[]{
                    new ElementRule(IntervalList.class, 1, Integer.MAX_VALUE)
                }),
                new ElementRule(POPULATION_TREE, new XMLSyntaxRule[]{
                    new ElementRule(TreeModel.class, 1, Integer.MAX_VALUE)
                })
            ),
            new ElementRule(GROUP_SIZES, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }, true),
            new ElementRule(SINGLE_BETA, new XMLSyntaxRule[] {
                    new ElementRule(Parameter.class),
            }, true),
            AttributeRule.newBooleanRule(RESCALE_BY_ROOT_ISSUE, true),
            AttributeRule.newBooleanRule(RANDOMIZE_TREE, true),
            AttributeRule.newBooleanRule(TIME_AWARE_SMOOTHING, true),
            AttributeRule.newBooleanRule(OLD_SKYRIDE, true),
            AttributeRule.newBooleanRule(BUILD_MAPPING, true)
    };

}